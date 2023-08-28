import argparse
import datetime as dt
import sys

import numpy as np
import xarray as xr
from dask.diagnostics import ProgressBar

from .common import encoding
from .rwguide.xarray import pvgradient, zonalization


def floats(arg):
    return np.asarray([float(x) for x in arg.split(",")])

parser = argparse.ArgumentParser()
parser.add_argument("--isentropes", type=floats, default=[330.], metavar="LVLS",
        help="isentropic levels to compute on (in K)")
parser.add_argument("--save-pv", action="store_true", help="output original PV fields")
parser.add_argument("--save-grad", action="store_true", help="output log(grad(PV)) fields")
parser.add_argument("--grad-exclude", type=float, default=0.1, metavar="PV",
        help="PV exclusion threshold for gradient computation (in PVU)")
parser.add_argument("--scale", type=float, default=120.,
        help="window width for rolling zonalization (in deg or km)")
parser.add_argument("--taper", type=float, default=0.8, help="alpha of Tukey window")
parser.add_argument("--fixed-km-scale", action="store_true",
        help="set if window width (--scale) is given in km instead of degrees longitude")
parser.add_argument("infiles", type=str, nargs="+")
parser.add_argument("outfile", type=str)

if __name__ == "__main__":
    args = parser.parse_args()

    uvt = xr.open_mfdataset(args.infiles, chunks={ "time": 8 })
    uvt = uvt.assign_coords({ "isentrope": args.isentropes })

    out = [] # output fields
    enc = {} # encoding for output fields

    isen = pvgradient.isob_to_isen_all(uvt)
    # Add original PV field to output if specified
    if args.save_pv:
        out.append(isen["pv"])
        enc["pv"] = encoding.preset.pv

    # Rolling-zonalized field in output
    lon = isen.coords["longitude"]
    lat = isen.coords["latitude"]
    if args.fixed_km_scale:
        win = zonalization.fixed_km_window(lon, lat, width=args.scale, window=("tukey", args.taper))
    else:
        win = zonalization.fixed_deg_window(lon, lat, width=args.scale, window=("tukey", args.taper))
    bgsg = zonalization.rolling_mean_background(win, isen["sg"])
    out.append(zonalization.zonalize_rolling(win, isen["av"], isen["sg"], bgsg))
    enc["pv_rz"] = encoding.preset.pv

    # Add grad-log waveguide metric to output if specified
    if args.save_grad:
        out.append(pvgradient.norm_grad_log_abs(pvrz, threshold=args.grad_exclude))
        enc["ngl_pv_rz"] = encoding.preset.pv_grad_log

    out = xr.merge(out)
    out.attrs = {
        "created": dt.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC"),
        "command": " ".join(sys.argv)
    }
    with ProgressBar():
        out.to_netcdf(args.outfile, encoding=enc)

