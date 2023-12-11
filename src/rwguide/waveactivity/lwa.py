import numpy as np
from ..bindings import *
from . import _ext


def local_wave_activity(lat, av, sg, pv_bg):
    lat, (nlat,), __ = require(lat, shape=(None,))
    av, (nlev, __, nlon), __ = require(av, shape=(None, nlat, None), fold=True)
    sg, (__  , __, __  ), __ = require(sg, shape=(nlev, nlat, nlon), fold=True)
    pv_bg, (__, __, __), out_shape = require(pv_bg, shape=(nlev, nlat, nlon), fold=True)
    out = empty(nlev, nlat, nlon)
    _ext.lib.local_wave_activity(
        as_pointer(lat),
        as_pointer(av),
        as_pointer(sg),
        as_pointer(pv_bg),
        nlev,
        nlat,
        nlon,
        as_pointer(out)
    )
    return out.reshape(out_shape)

