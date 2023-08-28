import numpy as np
from ..bindings import *
from . import _ext


def horizontal_gradient(lat, field):
    # Latitude in degrees
    lat, (nlat,), __ = require(lat, shape=(None,))
    # 3D input field
    field, (nlev, __, nlon), out_shape = require(field, shape=(None, nlat, None), fold=True)
    # Output fields: one per direction
    grad_lon = empty(nlev, nlat, nlon)
    grad_lat = empty(nlev, nlat, nlon)
    # Delegate computation to the C-extension
    status = _ext.lib.horizontal_gradient(
        as_pointer(lat),
        as_pointer(field),
        nlev,
        nlat,
        nlon,
        as_pointer(grad_lon),
        as_pointer(grad_lat)
    )
    assert status == 0, f"horizontal_gradient exited with status {status}"
    # Reshape output to match input temperature field
    return grad_lon.reshape(out_shape), grad_lat.reshape(out_shape)


def curl(lat, u, v):
    # Latitude in degrees
    lat, (nlat,), __ = require(lat, shape=(None,))
    # 3D input vector field (horizontal components)
    u, (nlev, __, nlon), out_shape = require(u, shape=(None, nlat, None), fold=True)
    v, (__  , __, __  ), __        = require(v, shape=(nlev, nlat, nlon), fold=True)
    # Output field
    curl = empty(nlev, nlat, nlon)
    # Delegate computation to the C-extension
    status = _ext.lib.curl(
        as_pointer(lat),
        as_pointer(u),
        as_pointer(v),
        nlev,
        nlat,
        nlon,
        as_pointer(curl),
    )
    assert status == 0, f"curl exited with status {status}"
    # Reshape output to match input temperature field
    return curl.reshape(out_shape)


def norm_grad_log_abs(lat, pv, threshold=0.):
    lat, (nlat,), __ = require(lat, shape=(None,))
    pv, (nlev, __, nlon), out_shape = require(pv, shape=(None, nlat, None), fold=True)
    out = empty(nlev, nlat, nlon)
    # Delegate computation ot the C-extension
    _ext.lib.norm_grad_log_abs(
        as_pointer(lat),
        as_pointer(pv),
        threshold,
        nlev,
        nlat,
        nlon,
        as_pointer(out)
    );
    return out.reshape(out_shape)

