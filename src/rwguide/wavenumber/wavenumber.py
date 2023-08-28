import numpy as np
from ..bindings import *
from . import _ext


def stationary_wavenumber(lat, u):
    lat, (nlat,), __ = require(lat, shape=(None,))
    u, (nlev, nlat, nlon), out_shape = require(u, shape=(None, nlat, None), fold=True)
    ks = np.empty_like(u)
    _ext.lib.stationary_wavenumber(
        as_pointer(lat),
        as_pointer(u),
        nlev,
        nlat,
        nlon,
        as_pointer(ks)
    )
    return ks.reshape(out_shape)


def stationary_wavenumber_squared(lat, u):
    lat, (nlat,), __ = require(lat, shape=(None,))
    u, (nlev, nlat, nlon), out_shape = require(u, shape=(None, nlat, None), fold=True)
    ks2 = np.empty_like(u)
    _ext.lib.stationary_wavenumber_squared(
        as_pointer(lat),
        as_pointer(u),
        nlev,
        nlat,
        nlon,
        as_pointer(ks2)
    )
    return ks2.reshape(out_shape)
