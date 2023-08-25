import numpy as np
import scipy.signal as sps

from ..bindings import *
from . import _ext


def require_kernel(x, **kwargs):
    if x is None:
        return require([1], **kwargs)
    elif isinstance(x, tuple) and len(x) >= 2 and isinstance(x[0], str):
        return require(getattr(sps, x[0])(*x[1:]), **kwargs)
    else:
        return require(x, **kwargs)

def mean_along_contour(xvalue, field_x, field_y, area_weights, lat_kernel=None, lon_kernel=None):
    field_x, (nlev, nlat, nlon), x_shape = require(field_x, shape=(None, None, None), fold=True)
    field_y, (__  , __  , __  ), __      = require(field_y, shape=(nlev, nlat, nlon), fold=True)
    area_weights, __, __ = require(area_weights, shape=(nlat,))
    lat_kernel, (nlat_kernel,), __ = require_kernel(lat_kernel, shape=(None,))
    lon_kernel, (nlon_kernel,), __ = require_kernel(lon_kernel, shape=(None,))
    out = empty(nlev, nlon)
    # Delegate computation to the C-extension
    _ext.lib.mean_along_contour(
        xvalue,
        as_pointer(field_x),
        as_pointer(field_y),
        as_pointer(area_weights),
        as_pointer(lat_kernel),
        as_pointer(lon_kernel),
        nlev,
        nlat,
        nlon,
        nlat_kernel,
        nlon_kernel,
        as_pointer(out)
    )
    out_shape = (*x_shape[:-2], nlon) # aggregation reduces along lat dimension
    return out.reshape(out_shape)

