import numpy as np
import xarray as xr 

from .common import get_names
from .. import wavenumber as _wavenumber


def stationary_wavenumber(da_u, names=None):
    lat, lon = get_names(names, "lat", "lon")
    return xr.apply_ufunc(
        _wavenumber.stationary_wavenumber,
        da_u[lat],
        da_u,
        input_core_dims=[
            [lat],
            [lat, lon]
        ],
        output_core_dims=[
            [lat, lon]
        ],
        vectorize=False,
        dask="parallelized",
        output_dtypes=[
            np.float64
        ]
    ).rename("ks")


def stationary_wavenumber_squared(da_u, names=None):
    lat, lon = get_names(names, "lat", "lon")
    return xr.apply_ufunc(
        _wavenumber.stationary_wavenumber_squared,
        da_u[lat],
        da_u,
        input_core_dims=[
            [lat],
            [lat, lon]
        ],
        output_core_dims=[
            [lat, lon]
        ],
        vectorize=False,
        dask="parallelized",
        output_dtypes=[
            np.float64
        ]
    ).rename("ks2")

