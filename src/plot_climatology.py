import argparse

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as art3d
from cartopy.util import add_cyclic_point
import shapely.geometry
from cartopy.feature import OCEAN, COASTLINE

from .common.seasons import select_season
from .common.plotting import (
    PLATE, add_box, add_map, finalize_map, cf_style_freq, feature_to_patches, draw_compass_3d
)


def plot_3d_waveguide(ax, freq, isens):
    # Background cross section
    cf = ax.contourf(
        freq.sel({ "longitude": -180 }).values,
        freq.coords["latitude"].values,
        freq.coords["isentrope"].values,
        zdir="x",
        offset=-180,
        levels=cf_style_freq["levels"],
        colors=cf_style_freq["colors"][1:],
        extend="max",
        zorder=10
    )
    # Filled contours for each level
    for isen in isens:
        ax.contourf(
            freq.coords["longitude"].values,
            freq.coords["latitude"].values,
            freq.sel({ "isentrope": isen }).values,
            zdir="z",
            offset=isen,
            levels=cf_style_freq["levels"][1:],
            colors=cf_style_freq["colors"][2:],
            extend="max",
            zorder=isen
        )
    # Highlight and shadow for 320K isentrope
    for offset in [291, 320]:
        ax.contour(
            freq.coords["longitude"].values,
            freq.coords["latitude"].values,
            freq.sel({ "isentrope": 320 }).values,
            zdir="z",
            offset=offset,
            levels=[30, 50, 70],
            colors="black",
            linestyles="solid",
            linewidths=0.9,
            zorder=offset
        )
    # Highlight and shadow for 340K isentrope
    for offset in [291, 340]:
        ax.contour(
            freq.coords["longitude"].values,
            freq.coords["latitude"].values,
            freq.sel({ "isentrope": 340 }).values,
            zdir="z",
            offset=offset,
            levels=[30, 50, 70],
            colors="black",
            linestyles="dashed",
            linewidths=0.9,
            zorder=offset
        )
    return cf

def ints(arg):
    return np.asarray([int(x) for x in arg.split(",")])

parser = argparse.ArgumentParser()
parser.add_argument("--scale", type=int, default=60, metavar="PATH")
parser.add_argument("--levels", type=ints, default=None, metavar="LVLS",
        help="isentropes to include in the structure plots")
parser.add_argument("--scale-cmp", type=int, default=90, metavar="PATH")
parser.add_argument("--level-cmp", type=int, default=330, metavar="LVL", help="isentropic level (in K)")
parser.add_argument("--season", choices=["ALL", "DJF", "MAM", "JJA", "SON", "winter", "summer"], default="DJF")
parser.add_argument("--threshold", type=float, default=1.0e-6)
parser.add_argument("outfile", type=str, nargs="?", default=None)

if __name__ == "__main__":
    args = parser.parse_args()

    fig = plt.figure(figsize=(10, 6.2))
    # Zonal spectrum
    zx = add_box(fig, (0.01, 0.77, 0.22, 0.19), label="a", title=f"60-30°N PV $k$-spectrum")
    # Two small maps in the first row to compare window widths
    ax1 = add_map(fig, (0.32, 0.74, 0.27, 0.22), label="b", title=f"rolling zonalized | {args.level_cmp} K | {args.scale}°")
    ax2 = add_map(fig, (0.63, 0.74, 0.27, 0.23), label="c", title=f"rolling zonalized | {args.level_cmp} K | {args.scale_cmp}°")
    # Vertical structure 3D-plots in a second row
    sx = add_box(fig, (0.01, -0.01, 0.38, 0.66), label="d", title=f"Southern Hemisphere: rolling zonalized | {args.scale}°",
            projection="3d", computed_zorder=False, anchor="NW")
    nx = add_box(fig, (0.50, -0.01, 0.38, 0.66), label="e", title=f"Northern Hemisphere: rolling zonalized | {args.scale}°",
            projection="3d", computed_zorder=False, anchor="NW")
    # Colorbar axes
    cx = fig.add_axes((0.93, 0.05, 0.01, 0.90))

    # 3D waveguide occurrence structure
    occ = xr.open_mfdataset([f"data/PVrz-{lvl}K-{args.scale}deg-occur.nc" for lvl in args.levels])
    occ = occ["pv_rz_occurrence"]
    occ = select_season(occ, args.season)
    occ = occ.sel({ "threshold": args.threshold })

    # 3D structure: Northern Hemisphere
    cf = plot_3d_waveguide(nx, occ.sel({ "latitude": slice(85, 5) }), args.levels)
    #cb = fig.colorbar(cf, cax=cx, orientation="horizontal", spacing="proportional")
    #cb.set_label("waveguide occurrence frequency [%]")
    draw_compass_3d(nx, 290, 80, 390, scale=1.0)
    # Take ocean feature from cartopy, cut out NH and add at 291K
    nh_poly = shapely.geometry.Polygon([[-180,  5], [180,  5], [180, 85], [-180, 85]])
    for patch in feature_to_patches(OCEAN, intersect_with=nh_poly, facecolor="#069", edgecolor="#000", alpha=0.4):
        nx.add_patch(patch)
        art3d.pathpatch_2d_to_3d(patch, z=291)
    nx.set_yticks([15, 30, 45, 60, 75])
    nx.set_yticklabels(["15°N", "", "45°N", "", "75°N"])
    nx.set_ylim(85, 5)
    nx.view_init(elev=25, azim=220, roll=0)

    # 3D structure: Southern Hemisphere
    plot_3d_waveguide(sx, occ.sel({ "latitude": slice(-5, -85) }), args.levels)
    draw_compass_3d(sx, 290, -80, 390, scale=1.0)
    # Take ocean feature from cartopy, cut out NH and add at 291K
    sh_poly = shapely.geometry.Polygon([[-180,  -5], [180,  -5], [180, -85], [-180, -85]])
    for patch in feature_to_patches(OCEAN, intersect_with=sh_poly, facecolor="#069", edgecolor="#000", alpha=0.4):
        sx.add_patch(patch)
        art3d.pathpatch_2d_to_3d(patch, z=291)
    sx.set_yticks([-15, -30, -45, -60, -75])
    sx.set_yticklabels(["15°S", "", "45°S", "", "75°S"])
    sx.set_ylim(-5, -85)
    sx.view_init(elev=25, azim=140, roll=0)

    for ax in [nx, sx]:
        ax.tick_params(axis="x", pad=1)
        ax.tick_params(axis="y", pad=-3)
        ax.set_xticks([-150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150])
        ax.set_xticklabels(["150°W", "", "90°W", "", "30°W", "", "30°E", "", "90°E", "", "150°E"])
        ax.set_zticks(args.levels)
        ax.set_xlim(180, -180)
        ax.set_zlim(291, 346)
        ax.set_box_aspect((2.0, 1.0, 1.1), zoom=1.2)
        ax.set_zmargin(0)

    nx.set_zticklabels([])
    sx.set_zticklabels([])
    fig.text(0.452, 0.424, "345 K", ha="center")
    fig.text(0.452, 0.396, "340 K", ha="center")
    fig.text(0.452, 0.368, "335 K", ha="center")
    fig.text(0.452, 0.340, "330 K", ha="center")
    fig.text(0.452, 0.312, "325 K", ha="center")
    fig.text(0.452, 0.284, "320 K", ha="center")
    fig.text(0.452, 0.256, "315 K", ha="center")

    # Waveguide occurrence window width comparison
    for ax, scale in zip([ax1, ax2], [args.scale, args.scale_cmp]):
        data = xr.Dataset({
            "occ": xr.open_dataarray(f"data/PVrz-{args.level_cmp}K-{scale}deg-occur.nc"),
            "avg": xr.open_dataarray(f"data/PVrz-{args.level_cmp}K-{scale}deg-mean.nc")
        })
        data = data.sel({ "isentrope": args.level_cmp })
        data = select_season(data, args.season)

        zz = data["occ"].sel({ "threshold": args.threshold })
        zz, xlon = add_cyclic_point(zz, data.coords["longitude"])
        ylat = data.coords["latitude"]
        cf = ax.contourf(xlon, ylat, zz, **cf_style_freq, transform=PLATE)
        zz = add_cyclic_point(data["avg"])
        ct = ax.contour(xlon, ylat, zz, transform=PLATE, colors="k", linewidths=0.6,
                linestyles="solid", levels=[-4.5, -2.5, -0.5, 0.5, 2.5, 4.5])

        finalize_map(ax)
        ax.plot([-180, 180], [0, 0], color="black", linewidth=1.5, transform=PLATE)

    cb = fig.colorbar(cf, cax=cx, orientation="vertical", spacing="proportional")
    cb.set_label("waveguide occurrence frequency [%]")

    # Zonal wavenumber spectra
    zspec = xr.Dataset({
        "rz60": xr.open_dataarray(f"data/PVrz-{args.level_cmp}K-{args.scale}deg-sprop.nc"),
        "rz90": xr.open_dataarray(f"data/PVrz-{args.level_cmp}K-{args.scale_cmp}deg-sprop.nc"),
        "rm14": xr.open_dataarray(f"data/PVrm-{args.level_cmp}K-sprop.nc"),
        "orig": xr.open_dataarray(f"data/PV-{args.level_cmp}K-sprop.nc"),
        "clim": xr.open_dataset("data/mean-isen-sprop.nc")["pv_zspec"]
    })
    zspec = zspec.sel({ "isentrope": args.level_cmp })
    zspec = select_season(zspec, args.season)
    zspec = zspec.sel({ "latitude": slice(60, 30) }).mean(dim="latitude") # TODO
    zspec = zspec * 1e-3

    k = zspec.coords["wavenumber"]
    zx.plot(k, zspec["orig"], color="#000", linewidth=2.0, label="inst.")
    zx.plot(k, zspec["clim"], color="#999", linewidth=2.0, label="clim.")
    zx.plot(k, zspec["rm14"], color="#369", linewidth=2.0, label="14d")
    zx.plot(k, zspec["rz60"], color="#F00", linewidth=2.0, label=f"{args.scale}°")
    zx.plot(k, zspec["rz90"], color="#900", linewidth=2.0, label=f"{args.scale_cmp}°")
    zx.legend(loc="upper right", fontsize="small", framealpha=1)
    zx.set_xticks([1, 2, 3, 4, 5, 6, 7])
    zx.set_xlim(1, 7.5)
    zx.set_ylim(0, 11)
    zx.set_ylabel("Power [$\mathrm{PVU}^2$]")
    zx.yaxis.tick_right()
    zx.yaxis.set_label_position("right")

    if args.outfile is None:
        plt.show()
    else:
        fig.savefig(args.outfile, dpi=150)

