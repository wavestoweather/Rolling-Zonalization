import argparse
import datetime as dt
import sys

import numpy as np
import xarray as xr
from dask.diagnostics import ProgressBar

from .common import encoding
from .rwguide.xarray import pvgradient


def floats(arg):
    return np.asarray([float(x) for x in arg.split(",")])

parser = argparse.ArgumentParser()
parser.add_argument("--isentropes", type=floats, default=[330.], metavar="LVLS",
        help="isentropic levels to compute on (in K)")
parser.add_argument("infiles", type=str, nargs="+")
parser.add_argument("outfile", type=str)

if __name__ == "__main__":
    args = parser.parse_args()

    uvt = xr.open_mfdataset(args.infiles, chunks={ "time": 8 })
    uvt = uvt.assign_coords({ "isentrope": args.isentropes })

    isen = pvgradient.isob_to_isen_all(uvt)

    out = isen[["pv"]]
    out.attrs = {
        "created": dt.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC"),
        "command": " ".join(sys.argv)
    }
    with ProgressBar():
        out.to_netcdf(args.outfile, encoding={ "pv": encoding.preset.pv })

