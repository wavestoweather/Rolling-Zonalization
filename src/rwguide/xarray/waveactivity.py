import numpy as np
import xarray as xr

from . import common as _common
from .. import waveactivity as _waveactivity


def local_wave_activity(da_av, da_sg, da_pv_bg=None, *, vectorize=True, names=None):
    """Finite-amplitude local wave activity (FALWA or just LWA).

    Generalized for localized zonalized background states.

    .. _Ghinassi et al. (2018): https://doi.org/10.1175/MWR-D-18-0068.1

    Parameters
    ----------
    da_av : xarray.DataArray
        Absolute vorticity in 1 / s.
    da_sg : xarray.DataArray
        Isentropic density in kg / K / m².
    da_pv_bg : xarray.DataArray, optional
        Background-state potential vorticity in PVU. If no background is
        specified, a hemispherically zonalized background state is computed
        from the input, resulting in local wave activity as defined by
        `Ghinassi et al. (2018)`_.
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
    # Fall back to Ghinassi et al. (2018) LWA if no background specified
    if da_pv_bg is None:
        from .zonalization import zonalize
        da_pv_bg = zonalize(
            np.cos(np.deg2rad(da_av.coords[lat])), # area weights for the sphere
            da_av,
            da_sg,
            da_sg.mean(dim=lon), # zonal-mean background isentropic density
            vectorize=vectorize,
            names=names
        )
    # Turn background state PV into lat-lon fields if only meridional profile
    # are given (e.g. from hemispheric zonalization)
    da_pv_bg = _common.require_lon(da_pv_bg, da_av.coords[lon])
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

