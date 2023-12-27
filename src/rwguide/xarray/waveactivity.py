import numpy as np
import xarray as xr

from . import common as _common
from .. import waveactivity as _waveactivity


def local_wave_activity(da_av, da_sg, da_pv_bg=None, *, vectorize=True, names=None):
    """Finite-amplitude local wave activity (FALWA or just LWA).

    LWA is computed with a boxcounting quadrature as::

                    a     ⎧Δϕ
        A(λ,ϕ) =  ──────  ⎪   (q(λ,ϕ+ϕ') − Q(λ,ϕ)) σ cos(ϕ+ϕ') dϕ'
                  cos(ϕ)  ⎭0

    Note how the background state PV Q(λ,ϕ) is allowed to depend on longitude,
    to allow for both hemispherically and localized zonalized PV as background
    states.

    .. note::
        The current implementation does not produce negative LWA when
        encountering ground-intersecting isentropes. The (Lagrangian)
        information about the location of the surface is also for the Eulerian
        parts of the integral, so that negative values cannot arise and LWA is
        strictly nonnegative. See footnote 3 of `Nakamura and Solomon (2011)`_.
        This should not affect the value of LWA as a Rossby wave amplitude
        diagnostic for the extratropical upper troposphere, as considered by
        `Ghinassi et al. (2018)`_.

    Barotropic local wave activity can be computed with this function by
    setting isentropic density to 1 everywhere (use
    :py:func:`xarray.ones_like` to generate `da_sg`).

    .. versionadded:: 1.2

    .. _Nakamura and Solomon (2011): https://doi.org/10.1175/2011JAS3685.1
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
        Local wave activity in m / s.
    """
    lat, lon, av, sg = _common.get_names(names, "lat", "lon", "av", "sg")
    # Fall back to Ghinassi et al. (2018) LWA if no background specified
    if da_pv_bg is None:
        from .zonalization import zonalize, area_weights
        da_pv_bg = zonalize(
            area_weights(da_av.coords[lat]), # area weights for the sphere
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

