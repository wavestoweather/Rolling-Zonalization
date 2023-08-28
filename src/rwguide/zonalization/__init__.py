"""Computation of zonalized (equivalent-latitude) background states."""

import numpy as np

from . import weighting
from .density import rolling_mean_background
from .zonalize import zonalize, zonalize_rolling


def fixed_deg_window(lon, lat, width=60., window=("tukey", 0.8)):
    nlon = lon.size
    dlon = lon[1] - lon[0]
    nwin = weighting.oddify(width // dlon)
    win_x = weighting.get_window(window, nwin, fold=(nlon, "max"))
    win_y = weighting.area_weights(lat)
    return weighting.product(win_x, win_y)


def fixed_km_window(lon, lat, width=8000., window=("tukey", 0.8)):
    nlon = lon.size
    nwin = weighting.wavelength_to_width(width * 1000., lat, nlon) # km to m
    nwin = weighting.clip_width(nwin, nmax=10*nlon)
    nwin = weighting.oddify(nwin)
    win_x = weighting.get_window(window, nwin, fold=(nlon, "max"))
    win_y = weighting.area_weights(lat)
    return weighting.product(win_x, win_y)

