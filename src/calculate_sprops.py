import argparse 
import datetime as dt
import sys

import numpy as np
import pandas as pd
import xarray as xr
from dask.diagnostics import ProgressBar

from .common.xarray import group_count


def zonal_spectral_power(x, ks):
    spec = np.fft.rfft(x, axis=-1)[...,ks]
    spec = np.abs(spec)**2
    return spec.astype(np.float32)


parser = argparse.ArgumentParser()
parser.add_argument("--hemisphere", choices=["north", "south"], default=None)
parser.add_argument("infiles", type=str, nargs="+")
parser.add_argument("outfile", type=str)

if __name__ == "__main__":
    args = parser.parse_args()

    data = xr.open_mfdataset(args.infiles, chunks={ "time": 20 })
    # Hemisphere selection
    if args.hemisphere == "north":
        data = data.sel({ "latitude": slice(90, 0) })
    elif args.hemisphere == "south":
        data = data.sel({ "latitude": slice(0, -90) })

    # Wavenumbers kept from the zonal spectrum
    waven = np.arange(1, 13)
    data = data.assign_coords({
        "wavenumber": ("wavenumber", waven)
    })

    out = {}
    for var in data.data_vars:
        # Zonal spectral power
        zspec = xr.apply_ufunc(
            zonal_spectral_power,
            data[var],
            data.coords["wavenumber"],
            input_core_dims=[["longitude"], ["wavenumber"]],
            output_core_dims=[["wavenumber"]],
            vectorize=False, # use built-in vectorization of FFT
            dask="parallelized",
            output_dtypes=[np.float32]
        )

        # Save seasonal and all-time means if possible to compute, otherwise
        # just save the raw data
        if "time" in zspec.coords:
            by_season = zspec.groupby(zspec.coords["time"].dt.season)
            counts = group_count(by_season, "season")
            zspec = by_season.mean(dim="time")
            zspec_all = zspec.weighted(counts).mean(dim="season").assign_coords({ "season": "ALL" })
            out[f"{var}_zspec"] = xr.concat([zspec, zspec_all], dim="season")
        else:
            out[f"{var}_zspec"] = zspec

    out = xr.Dataset(out)
    out.attrs = {
        "created": dt.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC"),
        "command": " ".join(sys.argv)
    }
    with ProgressBar():
        out.to_netcdf(args.outfile)

