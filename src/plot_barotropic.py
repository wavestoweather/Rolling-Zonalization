import argparse

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from cartopy.util import add_cyclic_point

from .common.seasons import select_season
from .common.plotting import PLATE, add_box, add_map, finalize_map, cf_style_baro
from .waveguide.xarray import wavenumber, pvgradient

parser = argparse.ArgumentParser()
parser.add_argument("--level", type=float, default=300., metavar="LVL", help="vertical level in hPa or K")
parser.add_argument("--isen", action="store_true", help="set if data contains isentropic levels")
parser.add_argument("--season", choices=["ALL", "DJF", "MAM", "JJA", "SON", "winter", "summer"], default="DJF")
parser.add_argument("infile", type=str, default="data/ERA5-mean-isob.nc", metavar="PATH", help="path to mean data")
parser.add_argument("outfile", type=str, nargs="?", default=None)


if __name__ == "__main__":
    args = parser.parse_args()

    fig = plt.figure(figsize=(10, 2.4))

    level = "{:.0f} {}".format(args.level, ("K" if args.isen else "hPa"))

    ax1 = add_map(fig, (0.005, 0.31, 0.31, 0.675), label="a")
    ax2 = add_map(fig, (0.345, 0.31, 0.31, 0.675), label="b")
    ax3 = add_map(fig, (0.685, 0.31, 0.31, 0.675), label="c")

    cx1 = fig.add_axes((0.03, 0.195, 0.26, 0.051))
    cx2 = fig.add_axes((0.37, 0.195, 0.26, 0.051))
    cx3 = fig.add_axes((0.71, 0.195, 0.26, 0.051))

    # Update of the Hoskins and Ambrizzi (1993) plots
    clim = xr.open_dataset(args.infile)
    clim = clim.sel({ ("isentrope" if args.isen else "level"): args.level }, drop=True)
    cori = 2. * 7.292e-5 * np.sin(np.deg2rad(clim["latitude"]))
    avort = (cori + pvgradient.curl(clim["u"], clim["v"])).rename("av")
    ks = wavenumber.stationary_wavenumber(clim["u"])
    avort_grad = pvgradient.horizontal_gradient(avort)
    clim = xr.merge([clim["u"], avort, avort_grad, ks])
    clim = select_season(clim, args.season)

    lon = clim["longitude"].values
    lat = clim["latitude"].values

    # Zonal wind
    zz, lon = add_cyclic_point(clim["u"], lon)
    cf1 = ax1.contourf(
        lon, lat, zz,
        levels=[10, 15, 20, 25, 30, 35, 40, 45, 50],
        **cf_style_baro,
        transform=PLATE
    )
    cb = fig.colorbar(cf1, cax=cx1, orientation="horizontal", spacing="proportional")
    cb.set_label(r"zonal wind $u$ $[\mathrm{m\,s^{-1}}]$")

    # Gradient of absolute vorticity
    cf2 = ax2.contourf(
        lon, lat, add_cyclic_point(1e12 * np.sqrt(clim["gradx_av"]**2 + clim["grady_av"]**2)),
        levels=[20, 25, 30, 35, 40, 45, 50, 55, 60],
        **cf_style_baro,
        transform=PLATE
    )
    cb = fig.colorbar(cf2, cax=cx2, orientation="horizontal", spacing="proportional")
    cb.set_label(r"abs. vort. grad. $\| \nabla \zeta_a \|$ $[10^{-12}\;\mathrm{m^{-1} \, s^{-1}}]$")

    # Stationary wavenumber
    cf3 = ax3.contourf(
        lon, lat, add_cyclic_point(clim["ks"]),
        levels=[0, 3, 4, 5, 6, 7, 8, 10, 15],
        **cf_style_baro,
        transform=PLATE
    )
    cb = fig.colorbar(cf3, cax=cx3, orientation="horizontal", spacing="proportional")
    cb.set_label(r"stationary wavenumber $K_s$")

    # Final touches
    for ax in [ax1, ax2, ax3]:
        finalize_map(ax)
        ax.plot([-180, 180], [0, 0], color="black", linewidth=1.5, transform=PLATE)

    if args.outfile is None:
        plt.show()
    else:
        fig.savefig(args.outfile, dpi=150)

