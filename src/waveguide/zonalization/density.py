import numpy as np
from ..bindings import *
from . import _ext


def rolling_mean_background(window, field):
    window, (nlat, nwin), __ = require(window, shape=(None, None))
    field, (nlev, __, nlon), out_shape = require(field, shape=(None, nlat, None), fold=True)
    out = empty(nlev, nlat, nlon)
    # Delegate computation to the C-extension
    status = _ext.lib.rolling_mean_background(
        as_pointer(window),
        as_pointer(field),
        nlev,
        nlat,
        nlon,
        nwin,
        as_pointer(out)
    )
    assert status == 0, f"AHHHH"
    # Reshape output to match input temperature field
    return out.reshape(out_shape)


def zonal_mean_density(sg, weights=None):
    """Weighted zonal mean isentropic density"""
    if weights is None:
        return sg.mean(axis=-1)
    raise NotImplementedError()

