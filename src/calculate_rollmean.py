import argparse
import datetime as dt
import sys

import xarray as xr
from scipy.signal import get_window
from dask.diagnostics import ProgressBar

from .common import encoding


def rolling_mean(da, window, name):
    window = window / window.sum() # normalize
    window = xr.DataArray(window, dims=["window"])
    return da.rolling({ "time": window.size }, center=True).construct("window").dot(window).rename(name)


parser = argparse.ArgumentParser()
parser.add_argument("--variables", type=lambda arg: arg.split(","), default=None, help="variables to process")
parser.add_argument("--length", type=int, default=57, help="timesteps in windowg mean")
parser.add_argument("--window", type=str, default="boxcar", help="window type")
parser.add_argument("infiles", type=str, nargs="+")
parser.add_argument("outfile", type=str)

if __name__ == "__main__":
    args = parser.parse_args()

    data = xr.open_mfdataset(args.infiles, chunks={ "time": 8 })

    if args.variables is not None:
        data = data[args.variables]

    name = "{}_rm".format

    out = []
    # Compute the rolling temporal mean from the output PV field
    weights = get_window(args.window, args.length, fftbins=False)
    for var, da in data.data_vars.items():
        out.append(rolling_mean(da, weights, name(var)))

    enc = { name(var): enc_dct for var, enc_dct in encoding.extract_encoding(data).items() }

    out = xr.merge(out)
    out.attrs = {
        "created": dt.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC"),
        "command": " ".join(sys.argv)
    }
    with ProgressBar():
        out.to_netcdf(args.outfile, encoding=enc)

