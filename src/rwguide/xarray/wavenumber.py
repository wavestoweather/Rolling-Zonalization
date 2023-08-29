import numpy as np
import xarray as xr 

from .common import get_names
from .. import wavenumber as _wavenumber


def stationary_wavenumber(da_u, vectorize=True, names=None):
    """Barotropic stationary wavenumber.

    Computed as ``Kₛ = Re(sqrt(Kₛ²)) - Im(sqrt(Kₛ²))``, i.e. positive values of
    ``Kₛ²`` remain positive after taking the square root and negative values
    remain negative.

    See also: :py:func:`stationary_wavenumber_squared`.

    Parameters
    ----------
    da_u : xarray.DataArray
        Zonal wind component. Core dimensions: latitude, longitude.
    vectorize : bool, optional
        Use vectorization of :py:func:`xarray.apply_ufunc`.
    names : dict, optional
        Variable name override.

    Returns
    -------
    xarray.DataArray
        Stationary wavenumber.
    """
    lat, lon, ks = get_names(names, "lat", "lon", "ks")
    return xr.apply_ufunc(
        _wavenumber.stationary_wavenumber,
        da_u[lat],
        da_u,
        input_core_dims=[
            [lat],
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
    ).rename(ks)


def stationary_wavenumber_squared(da_u, vectorize=True, names=None):
    """Square of the barotropic stationary wavenumber.

    ``Kₛ² = a cos(ϕ)² u⁻¹ ∂ζₐ/∂ϕ``, where ``a = 6371.2 km`` is the radius of
    Earth, ``ϕ`` is latitude, ``u`` is zonal wind and ``ζₐ`` is absolute
    vorticity (computed from the zonal wind only). The implementation is based
    on a three-point finite difference stencil. ``Kₛ²`` is set to 0 at the
    northern and southern boundaries of the domain.

    See, e.g.: `Karoly (1983)`_, `Hoskins and Ambrizzi (1993)`_.

    .. _Karoly (1983): https://dx.doi.org/10.1016/0377-0265(83)90013-1
    .. _Hoskins and Ambrizzi (1993): https://dx.doi.org/10.1175/1520-0469(1993)050<1661:RWPOAR>2.0.CO;2

    Parameters
    ----------
    da_u : xarray.DataArray
        Zonal wind component in ``m/s``. Core dimensions: latitude, longitude.
    vectorize : bool, optional
        Use vectorization of :py:func:`xarray.apply_ufunc`.
    names : dict, optional
        Variable name override.

    Returns
    -------
    xarray.DataArray
        Stationary wavenumber squared.
    """
    lat, lon, ks = get_names(names, "lat", "lon", "ks")
    return xr.apply_ufunc(
        _wavenumber.stationary_wavenumber_squared,
        da_u[lat],
        da_u,
        input_core_dims=[
            [lat],
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
    ).rename(f"{ks}2")

