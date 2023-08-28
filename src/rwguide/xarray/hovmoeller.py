import numpy as np
import xarray as xr

from . import common as _common
from .. import hovmoeller as _hovmoeller


def mean_along_contour(xvalue, da_x, da_y, da_area=None, lat_kernel=None, lon_kernel=None, vectorize=True, names=None):
    lat, lon = _common.get_names(names, "lat", "lon")
    if da_area is None:
        da_area = np.cos(np.deg2rad(da_x.coords[lat]))
    return xr.apply_ufunc(
        _hovmoeller.mean_along_contour,
        xvalue,
        da_x,
        da_y,
        da_area,
        kwargs={
            "lat_kernel": lat_kernel,
            "lon_kernel": lon_kernel
        },
        input_core_dims=[
            [],
            [lat, lon],
            [lat, lon],
            [lat]
        ],
        output_core_dims=[
            [lon]
        ],
        vectorize=vectorize,
        dask="parallelized",
        output_dtypes=[
            np.float64
        ]
    ).rename(da_y.name)

