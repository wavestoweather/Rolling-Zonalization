import argparse
import datetime as dt
import sys

import numpy as np
import xarray as xr
import scipy.ndimage
from dask.diagnostics import ProgressBar

from .common.xarray import group_count
from .waveguide.xarray.pvgradient import norm_grad_log_abs


def mask_thresholds(ngl, thresholds):
    return ngl > thresholds[:,None,None]

def only_circumglobal(mask):
    labelled, n = scipy.ndimage.label(mask * 1)
    # Check which objects wrap the globe
    wrapping = set(labelled[:,0]) & set(labelled[:,-1])
    # Numpy just wraps the set in an object-array, need to convert to list
    wrapping = np.asarray(list(wrapping))
    # 0-label is the background
    wrapping = wrapping[wrapping > 0]
    # New mask with only the circumglobal objects
    return np.isin(labelled, wrapping)


def floats(arg):
    return np.asarray([float(a) for a in arg.split(",")])

DEFAULT_THRESHOLDS = [0.6e-6, 0.8e-6, 1.0e-6, 1.2e6, 1.4e-6]

parser = argparse.ArgumentParser()
parser.add_argument("--thresholds", type=floats, default=DEFAULT_THRESHOLDS, metavar="GRAD",
        help="waveguide detection thresholds for grad(log(PV)) field")
parser.add_argument("--grad-exclude", type=float, default=0.1, metavar="PV",
        help="PV exclusion threshold for gradient computation (in PVU)")
parser.add_argument("--filter", choices=["circumglobal"], default=None,
        help="Only count waveguides that fulfill the criterion")
parser.add_argument("infiles", type=str, nargs="+")
parser.add_argument("outfile", type=str)

if __name__ == "__main__":
    args = parser.parse_args()

    # Isentropes in descending order (top-down)
    infiles = list(reversed(sorted(args.infiles)))
    data = xr.open_mfdataset(
        infiles,
        combine="nested",
        concat_dim="isentrope",
        chunks={ "time": 40, "isentrope": 1 }
    )


    out = {}
    for var in data.data_vars:
        waveguide = norm_grad_log_abs(data[var], threshold=args.grad_exclude)
        wgn = waveguide.name
        # Basic waveguide detection with threshold
        waveguide = waveguide.to_dataset() # need a dataset to assign new dim
        waveguide = waveguide.assign_coords({ "threshold": args.thresholds })
        mask = xr.apply_ufunc(
            mask_thresholds,
            waveguide[wgn],
            waveguide.coords["threshold"],
            input_core_dims=[["latitude", "longitude"], ["threshold"]],
            output_core_dims=[["threshold", "latitude", "longitude"]],
            vectorize=True,
            dask="parallelized",
            output_dtypes=[bool]
        ).rename("occurrence")

        # Filtering
        if args.filter == "circumglobal":
            mask = xr.apply_ufunc(
                only_circumglobal,
                mask,
                input_core_dims=[["latitude", "longitude"]],
                output_core_dims=[["latitude", "longitude"]],
                vectorize=True,
                dask="parallelized",
                output_dtypes=[bool]
            )

        # Save seasonal and all-time means
        by_season = mask.groupby(waveguide["time"].dt.season)
        counts = group_count(by_season, "season")
        freq = 100. * by_season.mean(dim="time")
        freq_all = freq.weighted(counts).mean(dim="season").assign_coords({ "season": "ALL" })
        out[f"{var}_occurrence"] = xr.concat([freq, freq_all], dim="season")
        
    out = xr.Dataset(out)
    out.attrs = {
        "created": dt.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC"),
        "command": " ".join(sys.argv)
    }
    with ProgressBar():
        out.to_netcdf(args.outfile)

