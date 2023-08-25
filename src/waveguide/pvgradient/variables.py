import numpy as np
from ..bindings import *
from . import _ext


def _is_ascending(x):
    return np.all(np.diff(x) > 0)

def _is_descending(x):
    return np.all(np.diff(x) < 0)


def potential_temperature_isob(isob, t):
    # Pressure levels in hPa define the vertical coordinate
    isob, (nisob,), __ = require(isob, shape=(None,))
    # 3D temperature field in pressure coordinate
    t, (__, nlat, nlon), out_shape = require(t, shape=(nisob, None, None))
    # Output 3D potential temperature field
    th = empty(nisob, nlat, nlon)
    # Delegate computation to the C-extension
    status = _ext.lib.potential_temperature_isob(
        as_pointer(isob),
        as_pointer(t),
        nisob,
        nlat,
        nlon,
        as_pointer(th)
    )
    assert status == 0, f"potential_temperature_isob exited with status {status}"
    # Reshape output to match input temperature field
    return th.reshape(out_shape)


def isentropic_density_isob(isob, th):
    # Pressure levels in hPa define the vertical coordinate
    isob, (nisob,), __ = require(isob, shape=(None,))
    # 3D potential temperature field in pressure coordinate
    th, (__, nlat, nlon), out_shape = require(th, shape=(nisob, None, None))
    # Output 3D isentropic density field in pressure coordinate
    sg = empty(nisob, nlat, nlon)
    # Delegate to the C-extension
    status = _ext.lib.isentropic_density_isob(
        as_pointer(isob),
        as_pointer(th),
        nisob,
        nlat,
        nlon,
        as_pointer(sg)
    )
    assert status == 0, f"isentropic_density_isob exited with status {status}"
    # Reshape output to match input potential temperature field
    return sg.reshape(out_shape)


def mask_underground_isob(isob, psfc, field, fill=np.nan, inplace=True):
    # Pressure levels in hPa define the vertical coordinate
    isob, (nisob,), __ = require(isob, shape=(None,))
    # Surface pressure field defines boundary below which values are masked
    psfc, (nlat, nlon), __ = require(psfc, shape=(None, None))
    # Input field to be masked
    field, __, out_shape = require(field, shape=(nisob, nlat, nlon))
    # Work on a copy of the input field if masking not applied in-place
    if not inplace:
        field = field.copy()
    # Delegate to the C-extension
    status = _ext.lib.mask_underground_isob(
        as_pointer(isob),
        as_pointer(psfc),
        as_pointer(field),
        fill,
        nisob,
        nlat,
        nlon
    )
    assert status == 0, f"mask_underground_isob exited with status {status}"
    # Reshape masked field to match input field shape
    return field.reshape(out_shape)


def interpolate_isob_to_isen(isob, isen, th, *fields):
    # Pressure levels in hPa define the input vertical coordinate
    isob, (nisob,), __ = require(isob, shape=(None,))
    assert _is_ascending(isob), "pressure levels must be given top-down (ascending)"
    # Isentropic levels in K define the output vertical coordinate
    isen, (nisen,), __ = require(isen, shape=(None,))
    assert _is_descending(isen), "isentropic levels must be given top-down (descending)"
    # 3D potential temperature field in pressure coordinate
    th, (nisob, nlat, nlon), __ = require(th, shape=(nisob, None, None))
    # Interpolated output fields are handled as variadic arguments
    var_args = []
    # Output 3D interpolated pressure on isentropes (convenient to include in
    # list of variadic arguments here)
    var_args.append(empty(nisen, nlat, nlon))
    # For every additional field that is interpolated, include an input and
    # output field argument
    for field in fields:
        field, __, __ = require(field, shape=(nisob, nlat, nlon))
        var_args.append(field)
        var_args.append(empty(nisen, nlat, nlon))
    nvar = len(fields)
    # Delegate to the C-extension
    status = _ext.lib.interpolate_isob_to_isen(
        as_pointer(isob),
        as_pointer(isen),
        as_pointer(th),
        nisob,
        nisen,
        nlat,
        nlon,
        nvar,
        *map(as_pointer, var_args)
    )
    assert status == 0, f"interpolate_isob_to_isen exited with status {status}"
    # Return interpolated output fields, reshaped to match shape of isentropes
    out_shape = (nisen, nlat, nlon)
    out_args = tuple(arg.reshape(out_shape) for arg in var_args[::2])
    return out_args if nvar > 0 else out_args[0]


def mask_underground_isen(p, psfc, field, fill=np.nan, inplace=True):
    # Isentropic pressure field provides reference values for masking with psfc
    p, (nisen, nlat, nlon), __ = require(p, shape=(None, None, None))
    # Surface pressure field defines boundary below which values are masked
    psfc, __, __ = require(psfc, shape=(nlat, nlon))
    # Input field to be masked
    field, __, out_shape = require(field, shape=(nisen, nlat, nlon))
    # Work on a copy of the input field if masking not applied in-place
    if not inplace:
        field = field.copy()
    # Delegate to the C-extension
    status = _ext.lib.mask_underground_isen(
        as_pointer(p),
        as_pointer(psfc),
        as_pointer(field),
        fill,
        nisen,
        nlat,
        nlon
    )
    assert status == 0, f"mask_underground_isen exited with status {status}"
    # Reshape masked field to match input field shape
    return field.reshape(out_shape)


def absolute_vorticity(lat, u, v):
    # Latitude in degrees
    lat, (nlat,), __ = require(lat, shape=(None,))
    # Horizontal wind field
    u, (nlev, __, nlon), out_shape = require(u, shape=(None, nlat, None))
    v, __              , __        = require(v, shape=(nlev, nlat, nlon))
    av = empty(nlev, nlat, nlon)
    # Delegate computation to the C-extension
    status = _ext.lib.absolute_vorticity(
        as_pointer(lat),
        as_pointer(u),
        as_pointer(v),
        nlev,
        nlat,
        nlon,
        as_pointer(av)
    )
    assert status == 0, f"absolute_vorticity exited with status {status}"
    return av.reshape(out_shape)


def potential_vorticity_isen(av, sg):
    # In: isentropic density, absolute vorticity
    av, (nisen, nlat, nlon), out_shape = require(av, shape=(None , None, None))
    sg, __                 , __        = require(sg, shape=(nisen, nlat, nlon))
    # Out: potential vorticity
    pv = empty(nisen, nlat, nlon)
    # Delegate computation to the C-extension
    status = _ext.lib.potential_vorticity_isen(
        as_pointer(av),
        as_pointer(sg),
        nisen,
        nlat,
        nlon,
        as_pointer(pv)
    )
    return pv.reshape(out_shape)


def isob_to_isen_all(isob, isen, lat, t, u, v, psfc=None):
    # Pressure levels in hPa define the input vertical coordinate
    isob, (nisob,), __ = require(isob, shape=(None,))
    assert _is_ascending(isob), "pressure levels must be given top-down (ascending)"
    # Isentropic levels in K define the output vertical coordinate
    isen, (nisen,), __ = require(isen, shape=(None,))
    assert _is_descending(isen), "isentropic levels must be given top-down (descending)"
    # Latitude grid
    lat, (nlat,), __ = require(lat, shape=(None,))
    assert _is_descending(lat), "latitude grid must be given N-to-S (descending)"
    # In: u, v, t on pressure levels
    t, (nvec, __, __, nlon), in_shape = require(t, shape=(None, nisob, nlat, None), fold=True)
    v, __                  , __       = require(v, shape=(nvec, nisob, nlat, nlon), fold=True)
    u, __                  , __       = require(u, shape=(nvec, nisob, nlat, nlon), fold=True)
    # In: surface pressure field to mask underground regions
    if psfc is not None:
        psfc, __, __ = require(psfc, shape=(nvec, nlat, nlon))
    # Out: u, v, t, sg, av, pv on isentropic levels
    pp = empty(nvec, nisen, nlat, nlon)
    uu = empty(nvec, nisen, nlat, nlon)
    vv = empty(nvec, nisen, nlat, nlon)
    sg = empty(nvec, nisen, nlat, nlon)
    av = empty(nvec, nisen, nlat, nlon)
    pv = empty(nvec, nisen, nlat, nlon)
    # Delegate computation to the C-extension
    _ext.lib.isob_to_isen_all(
        as_pointer(isob),
        as_pointer(isen),
        as_pointer(lat),
        as_pointer(t),
        as_pointer(u),
        as_pointer(v),
        as_pointer(psfc),
        nvec,
        nisob,
        nisen,
        nlat,
        nlon,
        as_pointer(pp),
        as_pointer(uu),
        as_pointer(vv),
        as_pointer(sg),
        as_pointer(av),
        as_pointer(pv)
    )
    out_shape = in_shape[:-3] + (nisen, nlat, nlon)
    return (
        pp.reshape(out_shape),
        uu.reshape(out_shape),
        vv.reshape(out_shape),
        sg.reshape(out_shape),
        av.reshape(out_shape),
        pv.reshape(out_shape)
    )

