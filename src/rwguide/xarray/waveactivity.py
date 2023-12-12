import numpy as np
import xarray as xr

from . import common as _common
from .. import waveactivity as _waveactivity


def local_wave_activity(da_av, da_sg, da_pv_bg, *, vectorize=True, names=None):
    """Finite-amplitude local wave activity.

    Generalized for localized zonalized background states.

    Parameters
    ----------
    da_av : xarray.DataArray
        Absolute vorticity in 1 / s.
    da_sg : xarray.DataArray
        Isentropic density in kg / K / m².
    da_pv_bg : xarray.DataArray
        Background-state potential vorticity in PVU.
    vectorize : boolean, optional
        Use vectorization of :py:func:`xarray.apply_ufunc`.
    names : dict, optional
        Variable name override.

    Returns
    -------
    xarray.DataArray
        Local wave activity in m/s.
    """
    lat, lon, av, sg = _common.get_names(names, "lat", "lon", "av", "sg")
    # Turn background state PV into lat-lon fields if only meridional profile
    # are given (e.g. from hemispheric zonalization)
    if lon not in da_pv_bg.coords:
        da_pv_bg = da_pv_bg.expand_dims({ lon: da_av.coords[lon] }, axis=-1)
    return xr.apply_ufunc(
        _waveactivity.local_wave_activity,
        da_av[lat],
        da_av,
        da_sg,
        da_pv_bg,
        input_core_dims=[
            [lat],
            [lat, lon],
            [lat, lon],
            [lat, lon]
        ],
        output_core_dims=[
            [lat, lon]
        ],
        vectorize=vectorize,
        dask="parallelized",
        output_dtypes=[
            np.float64
        ]
    ).rename(f"lwa")

