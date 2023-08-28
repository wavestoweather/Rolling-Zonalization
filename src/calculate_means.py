import argparse
import datetime as dt
import sys

import numpy as np
import xarray as xr
from dask.diagnostics import ProgressBar

from .common.encoding import extract_encoding
from .common.xarray import group_count
from .waveguide.xarray.pvgradient import isob_to_isen_all


def floats(arg):
    return np.asarray([float(x) for x in arg.split(",")])

parser = argparse.ArgumentParser()
parser.add_argument("--levels", type=floats, default=None, metavar="LVLS",
        help="which levels to use (pressure in hPa or isentropes in K), omit to take from input data")
parser.add_argument("--isen", action="store_true", help="interpolate to isentropic levels")
parser.add_argument("infiles", type=str, nargs="+")
parser.add_argument("outfile", type=str)

if __name__ == "__main__":
    args = parser.parse_args()

    data = xr.open_mfdataset(args.infiles, chunks={ "time": 4 })
    levname = "level"
    encoding = extract_encoding(data)

    # Option to interpolate to isentropes
    if args.isen:
        assert args.levels is not None, "no isentropic levels to interpolate to specified"
        levname = "isentrope"
        data = data.assign_coords({ "isentrope": args.levels })
        data = isob_to_isen_all(data)
        encoding = None # TODO
    # Use specified pressure levels of the data only (all if none specified)
    elif args.levels is not None:
        data = data.sel({ levname: args.levels })

    by_season = data.groupby(data["time"].dt.season)
    counts = group_count(by_season, "season")
    # Compute mean of each season first, then combine into all-season
    average = by_season.mean(dim="time")
    average_all = average.weighted(counts).mean(dim="season").assign_coords({ "season": "ALL" })
    # Combine and add metadata
    out = xr.concat([average, average_all], dim="season")
    out.attrs = {
        "created": dt.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC"),
        "command": " ".join(sys.argv)
    }
    with ProgressBar():
        out.to_netcdf(args.outfile, encoding=encoding)

