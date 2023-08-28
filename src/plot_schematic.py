import argparse
import datetime as dt
import os

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from cartopy.util import add_cyclic_point

from .common.plotting import (
    PLATE, add_map, finalize_map, draw_arrow, cf_style_pv, cf_style_gl, cf_style_aw,
    ct_style_v, add_legend_v
)
from .waveguide.xarray import pvgradient, zonalization


def fill_window(window, lon):
    lat = window.coords["latitude"]
    nlon = window.coords["longitude"].size
    if nlon == lon.size:
        return window
    filled = np.zeros((lat.size, lon.size), dtype=window.dtype)
    l = (lon.size - nlon) // 2
    r = l + nlon
    filled[:,l:r] = window.values
    return xr.DataArray(filled, coords=[lat, lon], name=window.name, attrs=window.attrs)


parser = argparse.ArgumentParser()
parser.add_argument("--date", type=str, default="2016-09-26T00:00",
        help="analysis date to build the schematic from")
parser.add_argument("--scale", type=float, default=60.,
        help="window width for rolling zonalization (in deg or km)")
parser.add_argument("--taper", type=float, default=0.0, help="alpha of Tukey window")
parser.add_argument("--fixed-km-scale", action="store_true",
        help="set if window width (--scale) is given in km instead of degrees longitude")
parser.add_argument("--isentrope", type=float, default=330., metavar="LVL",
        help="isentropic level (in K)")
parser.add_argument("--grad-exclude", type=float, default=0.1, metavar="PV",
        help="PV exclusion threshold for gradient computation (in PVU)")
parser.add_argument("--mean-width", type=int, default=14, metavar="DAYS",
        help="window width for rolling mean comparison (in d)")
parser.add_argument("infile", type=str, default="data/ERA5/ERA5-2016-tuv-1.5.nc")
parser.add_argument("outfile", type=str, nargs="?", default=None)

if __name__ == "__main__":
    args = parser.parse_args()

    date = pd.to_datetime(args.date)
    data = xr.open_dataset(args.infile).sel({ "latitude": slice(90, 0) })
    # Add isentropic level for interpolation
    data = data.assign_coords({ "isentrope": [args.isentrope] })
    # Narrow down on selected date/interval
    date_delta = dt.timedelta(hours=args.mean_width*24) / 2
    interval = slice(date-date_delta, date+date_delta)
    data_mean = data.sel({ "time": interval })
    data_inst = data.sel({ "time": args.date }, drop=True)

    lat = data["latitude"]
    lon = data["longitude"]
    gridspacing = abs(lon.values[1] - lon.values[0])

    # Interpolate data to selected isentropic level
    isen_mean = pvgradient.isob_to_isen_all(data_mean).squeeze().mean(dim="time")
    isen_inst = pvgradient.isob_to_isen_all(data_inst).squeeze()

    # PV gradient diagnostics for two background states:
    # 1) Temporal mean
    mean_grad = pvgradient.norm_grad_log_abs(isen_mean["pv"], threshold=0.1)
    # 2) Rolling zonalization
    if args.fixed_km_scale:
        win = zonalization.fixed_km_window(lon, lat, width=args.scale, window=("tukey", args.taper))
    else:
        win = zonalization.fixed_deg_window(lon, lat, width=args.scale, window=("tukey", args.taper))
    bgsg = zonalization.rolling_mean_background(win, isen_inst["sg"])
    znld = zonalization.zonalize_rolling(win, isen_inst["av"], isen_inst["sg"], bgsg)
    znld_grad = pvgradient.norm_grad_log_abs(znld, threshold=0.1)

    fig = plt.figure(figsize=(10, 6.2))
    maps = dict()
    # a) Input PV field
    maps["input"] = add_map(fig, (0.01, 0.77, 0.45, 0.20), label="a",
            title=f"{date:%Y-%m-%d %H:%M} UTC | {args.isentrope:.0f} K")
    # b) Individual in-window zonalizations
    maps["zonal4"] = add_map(fig, (0.52, 0.86, 0.200, 0.120), label="b")
    maps["zonal3"] = add_map(fig, (0.55, 0.80, 0.200, 0.120))
    maps["zonal2"] = add_map(fig, (0.58, 0.74, 0.200, 0.120))
    maps["zonal1"] = add_map(fig, (0.61, 0.68, 0.200, 0.120))
    # c) 14-d time average PV field
    maps["mean"] = add_map(fig, (0.01, 0.42, 0.450, 0.20), label="c")
    # d) IRWZ output PV field
    maps["output"] = add_map(fig, (0.54, 0.42, 0.45, 0.20), label="d")
    # e) Mean-based background state log(pv) gradient
    maps["mean_pvg"] = add_map(fig, (0.01, 0.15, 0.45, 0.20), label="e")
    # f) IRWZ-based background state log(pv) gradient
    maps["irwz_pvg"] = add_map(fig, (0.54, 0.15, 0.45, 0.20), label="f")

    # Axes for colorbars
    cx_pv = fig.add_axes((0.22, 0.08, 0.35, 0.02))
    cx_gl = fig.add_axes((0.64, 0.08, 0.35, 0.02))

    # Connectors in left column
    draw_arrow(fig, (0.18, 0.76), (0.18, 0.62)) # zonal* to output
    draw_arrow(fig, (0.18, 0.41), (0.18, 0.35)) # output to irwz_pvg
    # Connectors in right column
    draw_arrow(fig, (0.46, 0.86), (0.51, 0.86)) # input to zonal*
    draw_arrow(fig, (0.71, 0.68), (0.71, 0.62)) # zonal* to output
    draw_arrow(fig, (0.71, 0.41), (0.71, 0.35)) # output to irwz_pvg
    # Arrow labels
    fig.text(0.80, 0.82, "rolling zonalization\n(60° window)", fontsize="medium")
    fig.text(0.20, 0.68, "temporal average\n(14-day centered)", fontsize="medium")
    # Dot-dot-dot connectors to indicate there are more window placements than
    # shown between the zonal* maps
    fig.text(0.54, 0.853, "...", fontsize="large", rotation=-45)
    fig.text(0.57, 0.793, "...", fontsize="large", rotation=-45)
    fig.text(0.60, 0.733, "...", fontsize="large", rotation=-45)

    # Input PV field
    zz, lonx = add_cyclic_point(isen_inst["pv"].values, lon.values)
    cf_pv = maps["input"].contourf(lonx, lat, zz, **cf_style_pv, transform=PLATE)
    cb_pv = fig.colorbar(cf_pv, cax=cx_pv, orientation="horizontal", spacing="proportional")
    cb_pv.set_ticks(cf_style_pv["levels"][1::2])
    cb_pv.set_label(r"potential vorticity $q$ [PVU]")
    # Meridional wind
    zz, lonx = add_cyclic_point(isen_inst["v"].values, lon.values)
    maps["input"].contour(lonx, lat, zz, **ct_style_v, transform=PLATE)
    add_legend_v(fig, loc="lower left", frameon=False, alignment="left")

    # 14-day-mean PV field
    zz, lonx = add_cyclic_point(isen_mean["pv"].values, lon.values)
    maps["mean"].contourf(lonx, lat, zz, **cf_style_pv, transform=PLATE)

    # 14-day-mean log(PV) gradient field
    zz, lonx = add_cyclic_point(mean_grad, lon)
    maps["mean_pvg"].contourf(lonx, lat, zz, **cf_style_gl, transform=PLATE)

    # Rolling zonalization example PV fields in window
    for roll, zonal in zip([-90, -50, -10, 30], ["zonal4", "zonal3", "zonal2", "zonal1"]):
        win_example = fill_window(win, lon)
        win_example = win_example.roll({ "longitude": int(roll / gridspacing) })
        zon_example = zonalization.zonalize(win_example, isen_inst["av"], isen_inst["sg"])
        # Zonalized PV profile as filled contours in window
        zz = np.repeat(zon_example.values, lon.size).reshape((lat.size, lon.size))
        zz = np.where(win_example > 0, zz, np.nan)
        maps[zonal].contourf(lon, lat, zz, **cf_style_pv, transform=PLATE)
        # 1.5 PVU contour of original PV in window
        zz = np.where(win_example > 0, isen_inst["pv"].values, np.nan)
        maps[zonal].contour(lon, lat, zz, levels=[1.5], linewidths=0.6, colors="black", transform=PLATE)
        # Window boundary and center line
        maps[zonal].contour(lon, lat, win_example, levels=[0], linewidths=1, colors="black", transform=PLATE)
        maps[zonal].plot([roll, roll], [0, 90], color="black", linewidth=1.5, transform=PLATE)

    # IWRZ output PV field
    zz, lonx = add_cyclic_point(znld, lon)
    maps["output"].contourf(lonx, lat, zz, **cf_style_pv, transform=PLATE)

    # IWRZ output grad ln PV field
    zz, lonx = add_cyclic_point(znld_grad, lon)
    cf_gl = maps["irwz_pvg"].contourf(lonx, lat, zz, **cf_style_gl, transform=PLATE)
    cb_gl = fig.colorbar(cf_gl, cax=cx_gl, orientation="horizontal", spacing="proportional")
    cb_gl.set_ticks([0.8e-6, 1.2e-6, 1.6e-6, 2.0e-6, 3.0e-6])
    cb_gl.set_ticklabels(["0.8", "1.2", "1.6", "2", "3"])
    cb_gl.set_label(r"waveguidability proxy $\| \nabla \log(q) \|$ [$10^{-6}$]")

    # 1.5 PVU contour of original PV for reference
    cx_pv.axvline(1.5, color="black", linewidth=1)
    for name in ["input", "mean", "mean_pvg", "output", "irwz_pvg"]:
        maps[name].contour(lon, lat, isen_inst["pv"], levels=[1.5], linewidths=0.6, colors="black", transform=PLATE)

    for ax in maps.values():
        finalize_map(ax, region="NH", aspect=1.1)

    if args.outfile is None:
        plt.show()
    else:
        fig.savefig(args.outfile, dpi=150)

