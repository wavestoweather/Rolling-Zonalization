import argparse
import functools

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point

from .common.plotting import PLATE, add_map, finalize_map, add_box, draw_arrow
from .common.plotting import cf_style_pv, cf_style_gl, ct_style_v, add_legend_v
from .waveguide.xarray import hovmoeller, pvgradient


def open_and_merge_data(year, isentrope=330, scale=60, threshold=0.1):
    # PV and meridional wind on isentrope
    tuvp = xr.open_dataset(f"data/ERA5/ERA5-{year}-tuv-1.5.nc", chunks={ "time": 4 })
    tuvp = tuvp.assign_coords({ "isentrope": [isentrope] })
    isen = pvgradient.isob_to_isen_all(tuvp)
    # Rolling zonalized PV and grad-log
    pvrz = xr.open_dataarray(f"data/PVrz-{isentrope:.0f}K-{scale:.0f}deg.nc", chunks={ "time": 4 })
    pvrz = pvrz.sel({ "time": isen.coords["time"].values })
    grad = pvgradient.norm_grad_log_abs(pvrz, threshold=threshold)
    return xr.merge([isen["v"], isen["pv"], pvrz, grad], join="exact", combine_attrs="drop")


def plot_episode(data, hx, dates, contour=1.5):
    xlon = data.coords["longitude"]
    ylat = data.coords["latitude"]
    ytime = data.coords["time"]
    hov_kernels = { "lat_kernel": ("boxcar", 7), "lon_kernel": ("boxcar", 7) }
    v_hov = hovmoeller.mean_along_contour(contour, data["pv_rz"], data["v"], **hov_kernels)
    g_hov = hovmoeller.mean_along_contour(contour, data["pv_rz"], data["ngl_pv_rz"], **hov_kernels)
    cf = hx.contourf(xlon, ytime, g_hov, **cf_style_gl)
    hx.contour(xlon, ytime, v_hov, **ct_style_v)
    hx.yaxis.tick_right()
    hx.set_xticks([-150, -90, -30, 30, 90, 150])
    hx.set_xticklabels(["150°W", "90°W", "30°W", "30°E", "90°E", "150°E"])
    hx.set_xticks([-180, -120, -60, 0, 60, 120, 180], minor=True)
    hx.set_xlim(-180, 180)
    hx.set_title(rf"Refined Hovmöller | $330\;\mathrm{{K}}$ | $q_\mathrm{{rz}} = {contour:.1f}\;\mathrm{{PVU}}$", fontsize="medium")
    for date, mx in dates:
        hx.axhline(np.datetime64(date), color="k", linewidth=1)
        if mx is not None:
            d = data.sel({ "time": date }, drop=True)
            mx.contourf(xlon, ylat, d["ngl_pv_rz"], transform=PLATE, **cf_style_gl)
            mx.contour(xlon, ylat, d["pv_rz"], transform=PLATE, levels=[contour], colors="k", linewidths=0.8)
            mx.contour(xlon, ylat, d["v"], transform=PLATE, **ct_style_v)
            mx.set_title(f"{date} UTC", fontsize="medium")
            finalize_map(mx, region="NH*", aspect=1.1)
    return cf


parser = argparse.ArgumentParser()
parser.add_argument("--scale", type=int, default=60, metavar="DEG")
parser.add_argument("--isentrope", type=float, default=330., metavar="LVL")
parser.add_argument("--grad-exclude", type=float, default=0.1, metavar="PVU",
        help="PV exclusion threshold for gradient computation")
parser.add_argument("outfile", type=str, nargs="?", default=None)

if __name__ == "__main__":
    args = parser.parse_args()

    # Apply provided configuration for data
    open_data = functools.partial(
        open_and_merge_data,
        isentrope=args.isentrope,
        scale=args.scale,
        threshold=args.grad_exclude
    )

    fig = plt.figure(figsize=(10, 6.2))

    mxA1 = add_map(fig, (0.01, 0.82, 0.40, 0.14), label="a", title="")
    mxA2 = add_map(fig, (0.01, 0.63, 0.40, 0.14), label="b", title="")
    hx1 = add_box(fig, (0.45, 0.65, 0.35, 0.31), label="c", title="")
    
    mxB1 = add_map(fig, (0.01, 0.39, 0.40, 0.14), label="d", title="")
    mxB2 = add_map(fig, (0.01, 0.20, 0.40, 0.14), label="e", title="")
    mxB3 = add_map(fig, (0.01, 0.01, 0.40, 0.14), label="f", title="")
    hx2 = add_box(fig, (0.45, 0.04, 0.35, 0.49), label="g", title="")

    cx = fig.add_axes((0.93, 0.15, 0.01, 0.80))

    # Episode 1: 2016-12-15 to 2016-12-21
    data = open_data(2016).squeeze()
    data = data.sel({
        "latitude": slice(90, 0),
        "time": slice("2016-12-14 00:00", "2016-12-20 00:00")
    })
    plot_episode(data.compute(), hx1, [
        ("2016-12-17 12:00", mxA1),
        ("2016-12-15 12:00", mxA2)
    ])
    hx1.set_yticks([np.datetime64(f"2016-12-{d:0>2}") for d in [14, 15, 16, 17, 18, 19, 20]])
    hx1.set_yticklabels([
        "14 Dec",
        "15 Dec",
        "16 Dec",
        "17 Dec",
        "18 Dec",
        "19 Dec",
        "20 Dec 2016",
    ])

    # Episode 2: 2018-01-02 to 2018-01-15
    data = open_and_merge_data(2018).squeeze()
    data = data.sel({
        "latitude": slice(90, 0),
        "time": slice("2018-01-04 00:00", "2018-01-12 00:00")
    })
    cf = plot_episode(data.compute(), hx2, [
        ("2018-01-09 00:00", mxB1),
        ("2018-01-07 00:00", mxB2),
        ("2018-01-05 00:00", mxB3)
    ], contour=1.0)
    hx2.set_yticks([np.datetime64(f"2018-01-{d:0>2}") for d in [4, 5, 6, 7, 8, 9, 10, 11, 12]])
    hx2.set_yticklabels([
        "4 Jan",
        "5 Jan",
        "6 Jan",
        "7 Jan",
        "8 Jan",
        "9 Jan",
        "10 Jan",
        "10 Jan",
        "12 Jan 2018",
    ])

    cb = plt.colorbar(cf, cax=cx, orientation="vertical", spacing="proportional")
    cb.set_ticks(cf_style_gl["levels"])
    cb.set_ticklabels(["{:.1f}".format(lvl * 1e6).rstrip(".0") for lvl in cf_style_gl["levels"]])
    cb.set_label(r"waveguidability proxy $\| \nabla \log(q_\mathrm{rz}) \|$ [$10^{-6}$]")

    add_legend_v(fig, loc="lower right", frameon=False, fontsize="small", title_fontsize="small")

    draw_arrow(fig, (0.445, 0.836), (0.42, 0.865))
    draw_arrow(fig, (0.445, 0.735), (0.42, 0.715))
    draw_arrow(fig, (0.445, 0.350), (0.42, 0.390))
    draw_arrow(fig, (0.445, 0.225), (0.42, 0.240))
    draw_arrow(fig, (0.445, 0.100), (0.42, 0.090))

    if args.outfile is None:
        plt.show()
    else:
        fig.savefig(args.outfile, dpi=150)

