import numpy as np
import xarray as xr

from . import common as _common
from .. import zonalization as _zonalization


def _window(window_fun, da_lon, da_lat, name, **kwargs):
    values = window_fun(da_lon.values, da_lat.values, **kwargs)
    nlon = values.shape[1]
    coords = {
        da_lat.name: da_lat.values,
        da_lon.name: da_lon.values[:nlon] - da_lon.values[0]
    }
    return xr.DataArray(values, coords=coords, name=name, attrs=kwargs)
    

def fixed_km_window(da_lon, da_lat, name="fixed_km_window", **kwargs):
    """Weighting window with a constant-distance width (in kilometers).

    Combines a latitude-dependent cosine area weighting and a longitudinally
    varying weighting window from :py:mod:`scipy.signal`.

    The constant width in terms of actual distance should ensure a consistent
    cut-off wavelength across the sphere. The window is returned folded (see
    :py:func:`rwguide.zonalization.weighting.fold_periodic`).

    .. note::
        We found in testing that a constant-distance width window has rather
        poor PV conservation properties when used with rolling zonalization.
        Consider using a constant-longitude window instead
        (:py:func:`fixed_deg_window`).

    Parameters
    ----------
    da_lon : xarray.DataArray
        Longitude coordinate.
    da_lat : xarray.DataArray
        Latitude coordinate.
    name : str, optional
        Name assigned to the output *DataArray*.
    kwargs : any
        Keyword arguments given to :py:func:`rwguide.zonalization.fixed_km_window`.

    Returns
    -------
    xarray.DataArray
        Weighting window.
    """
    return _window(_zonalization.fixed_km_window, da_lon, da_lat, name, **kwargs)


def fixed_deg_window(da_lon, da_lat, name="fixed_deg_window", **kwargs):
    """Weighting window with a constant-longitude width (in degrees longitude).

    Combines a latitude-dependent cosine area weighting and a longitudinally
    varying weighting window from :py:mod:`scipy.signal`.

    The constant width in terms of degrees longitude should ensure a consistent
    cut-off wavenumber across the sphere.

    Parameters
    ----------
    da_lon : xarray.DataArray
        Longitude coordinate.
    da_lat : xarray.DataArray
        Latitude coordinate.
    name : str, optional
        Name assigned to the output *DataArray*.
    kwargs : any
        Keyword arguments given to :py:func:`rwguide.zonalization.fixed_deg_window`.

    Returns
    -------
    xarray.DataArray
        Weighting window.
    """
    return _window(_zonalization.fixed_deg_window, da_lon, da_lat, name, **kwargs)


def rolling_mean_background(da_area, da_sg, *, vectorize=True, names=None):
    """Weighted longitudinal rolling mean.

    A simple procedure to approximate the background-state isentropic density
    field for rolling zonalization.

    Parameters
    ----------
    da_area : xarray.DataArray
        Weighting window.
    da_sg : xarray.DataArray
        Isentropic density field. Core dimensions: latitude, longitude.
    vectorize : boolean, optional
        Use vectorization of :py:func:`xarray.apply_ufunc`.
    names : dict, optional
        Variable name override.

    Returns
    -------
    xarray.DataArray
        Rolling-mean background isentropic density.
    """
    lat, lon, sg = _common.get_names(names, "lat", "lon", "sg")
    # To avoid a coordinate name conflict, rename the longitude dimension of
    # the area weights (they are relative)
    lon_area = f"{lon}_area"
    return xr.apply_ufunc(
        _zonalization.rolling_mean_background,
        da_area.rename({ lon: lon_area }),
        da_sg,
        input_core_dims=[
            [lat, lon_area],
            [lat, lon],
        ],
        output_core_dims=[
            [lat, lon]
        ],
        vectorize=vectorize,
        dask="parallelized",
        output_dtypes=[
            np.float64
        ]
    ).rename(f"{sg}_rm")


def zonalize(da_area, da_av, da_sg, da_bg=None, *, vectorize=True, names=None):
    """Zonalization (hemispheric or sectoral/fixed-window).

    Integrations are carried out with a conditional boxcounting
    quadrature scheme. Northern and southern hemispheres are automatically
    detected and zonalized separately. Regions with 0-valued isentropic density
    are omitted in the surface integrals.

    The PV contours for the zonalizations are automatically determined based on
    the input data: 10 contours, linearly spaced between 0 and the maximum (NH)
    or minimum (SH) value of PV found in the domain, are first zonalized to
    provide an approximation of the zonalized PV profile. The first guess PV
    profile is then interpolated to the input latitude grid and the
    interpolated PV values at each grid latitude are then used as contours in
    a second zonalization pass to obtain a refined zonalized profile. The
    refined zonalized PV profile is finally interpolated to the input latitude
    and returned.

    E.g., hemispheric zonalization based on a zonal-mean background state isentropic
    density profile `(Ghinassi et al. 2018)`_:

    >>> area = np.cos(np.deg2rad(isen.coords["latitude"]))
    >>> zonalize(area, isen["av"], isen["sg"], isen["sg"].mean(dim="longitude"))

    .. _(Ghinassi et al. 2018): https://doi.org/10.1175/MWR-D-18-0068.1

    Parameters
    ----------
    da_area : xarray.DataArray
        Weighting window or meridional profile of weighting coefficients (for
        "normal" hemispheric zonalization).
    da_av : xarray.DataArray
        Isentropic absolute vorticity field in 1 / s. Core dimensions:
        latitude, longitude.
    da_sg : xarray.DataArray
        Isentropic density field in kg / K / m². Core dimensions: latitude,
        longitude.
    da_bg : xarray.DataArray, optional
        Background state isentropic density field in kg / K / m². Core
        dimensions: latitude, longitude. If none given, *da_sg* is used for the
        background too. If no longitude dimension is found, the meridional
        profile is used at every longitude.
    vectorize : boolean, optional
        Use vectorization of :py:func:`xarray.apply_ufunc`.
    names : dict, optional
        Variable name override.

    Returns
    -------
    xarray.DataArray
        Zonalized PV in PVU.
    """
    lat, lon, pv = _common.get_names(names, "lat", "lon", "pv")
    if da_bg is None:
        da_bg = da_sg
    # Density background and area weights must be lat-lon fields
    da_bg = _common.require_lon(da_bg, da_av.coords[lon])
    da_area = _common.require_lon(da_area, da_av.coords[lon])
    # Determine indices where hemispheres can be separated
    nh_last = da_area[lat].sel({ lat: slice(90, 0) }).size
    sh_first = da_area[lat].size - da_area[lat].sel({ lat: slice(0, -90) }).size
    return xr.apply_ufunc(
        _zonalization.zonalize,
        da_area,
        da_av,
        da_sg,
        da_bg,
        kwargs={
            "nh_last": nh_last,
            "sh_first": sh_first
        },
        input_core_dims=[
            [lat, lon],
            [lat, lon],
            [lat, lon],
            [lat, lon],
        ],
        output_core_dims=[
            [lat]
        ],
        vectorize=vectorize,
        dask="parallelized",
        output_dtypes=[
            np.float64
        ]
    ).rename(f"{pv}_z")


def zonalize_rolling(da_area, da_av, da_sg, da_bg=None, *, vectorize=True, names=None):
    """Rolling zonalization.

    The PV contours for the zonalization are automatically determined based on
    the input data (see :py:func:`zonalize`). Integrations are carried out with
    a conditional boxcounting quadrature scheme. Northern and southern
    hemispheres are automatically detected and zonalized separately. Regions
    with 0-valued isentropic density are omitted in the surface integrals.

    Parameters
    ----------
    da_area : xarray.DataArray
        Weighting window.
    da_av : xarray.DataArray
        Isentropic absolute vorticity field in 1 / s. Core dimensions:
        latitude, longitude.
    da_sg : xarray.DataArray
        Isentropic density field in kg / K / m². Core dimensions: latitude,
        longitude.
    da_bg : xarray.DataArray, optional
        Background state isentropic density field in kg / K / m². Core
        dimensions: latitude, longitude. If none given, *da_sg* is used for the
        background too.
    vectorize : boolean, optional
        Use vectorization of :py:func:`xarray.apply_ufunc`.
    names : dict, optional
        Variable name override.

    Returns
    -------
    xarray.DataArray
        Rolling zonalized PV in PVU.
    """
    lat, lon, pv = _common.get_names(names, "lat", "lon", "pv")
    if da_bg is None:
        da_bg = da_sg
    # Determine indices where hemispheres can be separated
    nh_last = da_area[lat].sel({ lat: slice(90, 0) }).size
    sh_first = da_area[lat].size - da_area[lat].sel({ lat: slice(0, -90) }).size
    # To avoid a coordinate name conflict, rename the longitude dimension of
    # the area weights (they are relative)
    lon_area = f"{lon}_area"
    return xr.apply_ufunc(
        _zonalization.zonalize_rolling,
        da_area.rename({ lon: lon_area }),
        da_av,
        da_sg,
        da_bg,
        kwargs={
            "nh_last": nh_last,
            "sh_first": sh_first
        },
        input_core_dims=[
            [lat, lon_area],
            [lat, lon],
            [lat, lon],
            [lat, lon],
        ],
        output_core_dims=[
            [lat, lon]
        ],
        vectorize=vectorize,
        dask="parallelized",
        output_dtypes=[
            np.float64
        ]
    ).rename(f"{pv}_rz")

