import numpy as np
from ..bindings import *
from . import _ext


def zonalize(area, av, sg, bg, nh_last=None, sh_first=None):
    area, (nlat, nlon), __ = require(area, shape=(None, None))
    av, (nlev, __, __), av_shape = require(av, shape=(None, nlat, nlon), fold=True)
    sg, (__  , __, __), __       = require(sg, shape=(nlev, nlat, nlon), fold=True)
    bg, (__  , __, __), __       = require(bg, shape=(nlev, nlat, nlon), fold=True)
    out = empty(nlev, nlat)
    # By default, compute only northern hemisphere
    if nh_last is None:
        nh_last = nlat-1 if sh_first is None else sh_first
    if sh_first is None:
        sh_first = nlat-1 if nh_last is None else nh_last
    # Delegate computation to the C-extension
    _ext.lib.zonalize(
        as_pointer(area),
        as_pointer(av),
        as_pointer(sg),
        as_pointer(bg),
        nh_last,
        sh_first,
        nlev,
        nlat,
        nlon,
        as_pointer(out)
    )
    out_shape = av_shape[:-1] # zonalization reduces along lon dimension
    return out.reshape(out_shape)


def zonalize_rolling(area, av, sg, bg, nh_last=None, sh_first=None):
    area, (nlat, nwin), __ = require(area, shape=(None, None))
    av, (nlev, __, nlon), out_shape = require(av, shape=(None, nlat, None), fold=True)
    sg, (__  , __, __  ), __        = require(sg, shape=(nlev, nlat, nlon), fold=True)
    bg, (__  , __, __  ), __        = require(bg, shape=(nlev, nlat, nlon), fold=True)
    out = empty(nlev, nlat, nlon)
    # By default, compute only northern hemisphere
    if nh_last is None:
        nh_last = nlat-1 if sh_first is None else sh_first
    if sh_first is None:
        sh_first = nlat-1 if nh_last is None else nh_last
    # Delegate computation to the C-extension
    _ext.lib.zonalize_rolling(
        as_pointer(area),
        as_pointer(av),
        as_pointer(sg),
        as_pointer(bg),
        nh_last,
        sh_first,
        nlev,
        nlat,
        nlon,
        nwin,
        as_pointer(out)
    )
    return out.reshape(out_shape)

