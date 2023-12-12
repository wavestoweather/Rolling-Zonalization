import numpy as np
import xarray as xr

from .common import NAMES, get_names
from .. import pvgradient as _pvgradient


def horizontal_gradient(da, *, vectorize=True, names=None):
    """Horizontal gradient ∇· on the sphere.

    The sphere has the dimensions of Earth (6371.2 km radius).

    Parameters
    ----------
    da : xarray.DataArray
        Input field. Core dimensions: latitude, longitude.
    vectorize : boolean, optional
        Use vectorization of :py:func:`xarray.apply_ufunc`.
    names : dict, optional
        Variable name override.

    Returns
    -------
    xarray.Dataset
        Zonal (``gradx_*``) and meridional (``grady_*``) gradient of the scalar
        input field.
    """
    lat, lon = get_names(names, "lat", "lon")
    gx, gy = xr.apply_ufunc(
        _pvgradient.horizontal_gradient,
        da[lat],
        da,
        input_core_dims=[
            [lat],
            [lat, lon]
        ],
        output_core_dims=[
            [lat, lon],
            [lat, lon]
        ],
        vectorize=vectorize,
        dask="parallelized",
        output_dtypes=[
            np.float64,
            np.float64
        ]
    )
    gx = gx.rename(f"gradx_{da.name}")
    gy = gy.rename(f"grady_{da.name}")
    return xr.merge([gx, gy])


def curl(da_u, da_v, *, vectorize=True, names=None):
    """Vertical component of rotation ∇ ⨯ · on the sphere.

    The sphere has the dimensions of Earth (6371.2 km radius).

    Parameters
    ----------
    da_u : xarray.DataArray
        Input zonal vector component. Core dimensions: latitude, longitude.
    da_v : xarray.DataArray
        Input meridional vector component. Core dimensions: latitude,
        longitude.
    vectorize : boolean, optional
        Use vectorization of :py:func:`xarray.apply_ufunc`.
    names : dict, optional
        Variable name override.

    Returns
    -------
    xarray.DataArray
        Curl of the vector input field.
    """
    lat, lon = get_names(names, "lat", "lon")
    # TODO: pick a nice name based on the input names, e.g. u and v make relvort
    return xr.apply_ufunc(
        _pvgradient.curl,
        da_u[lat],
        da_u,
        da_v,
        input_core_dims=[
            [lat],
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
    ).rename(f"curl_{da_u.name}_{da_v.name}")


def norm_grad_log_abs(da, *, threshold=0., vectorize=True, names=None):
    """Magnitude of the gradient of the logarithm of the input field.

    Used as a waveguide diagnostic of isentropic Ertel PV. To first, order
    ∇(log(\|IPV\|)) ≈ 1/f₀ ∇QGPV ~ ∇²u, where IPV is isentropic potential
    vorticity, QGPV is quasi-geostropic potential vorticity and u is the
    zonal wind `(Martius et al. 2010)`_. 

    .. note::
        d/dx log(x/a) = d/dx log(x) = 1/x for all a = const, i.e. the result
        of this function does not depend on the unit of its argument.

    .. _(Martius et al. 2010): https://dx.doi.org/10.1175/2009JAS2995.1

    Parameters
    ----------
    da : xarray.DataArray
        Input field. Core dimensions: latitude, longitude.
    threshold : number, optional
        Threshold underneath which values of \|·\| are ignored, to deal with
        regions where the input field approaches zero and the logarithm tends
        toward negative infinity.
    vectorize : boolean, optional
        Use vectorization of :py:func:`xarray.apply_ufunc`.
    names : dict, optional
        Variable name override.

    Returns
    -------
    xarray.DataArray
        ‖∇(log(\|·\|))‖₂ of the input field in 1 / m.
    """
    lat, lon = get_names(names, "lat", "lon")
    return xr.apply_ufunc(
        _pvgradient.norm_grad_log_abs,
        da[lat],
        da,
        kwargs={
            "threshold": threshold
        },
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
    ).rename(f"ngl_{da.name}")


def potential_temperature_isob(da_t, *, names=None):
    """Potential temperature from temperature on pressure levels.

    The pressure coordinate (level) must be specified in hPa.

    Parameters
    ----------
    da_t : xarray.DataArray
        Temperature field in K. Core dimensions: level, latitude,
        longitude.
    names : dict, optional
        Variable name override.

    Returns
    -------
    xarray.DataArray
        Potential temperature in K.
    """
    isob, lat, lon, pt = get_names(names, "isob", "lat", "lon", "pt")
    return xr.apply_ufunc(
        _pvgradient.potential_temperature_isob,
        da_t[isob],
        da_t,
        input_core_dims=[
            [isob],
            [isob, lat, lon]
        ],
        output_core_dims=[
            [isob, lat, lon]
        ],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[
            np.float64
        ]
    ).rename(pt)


def isentropic_density_isob(da_pt, *, names=None):
    """Isentropic density from potential temperature on pressure levels.

    Isentropic density σ = - g⁻¹ ∂p/∂θ computed as -(g ∂θ/∂p)⁻¹ using
    a 2nd-order finite difference stencil. 0-values are assigned to grid points
    with unstable stratification.

    The pressure coordinate (level) must be specified in hPa.

    Parameters
    ----------
    da_pt : xarray.DataArray
        Potential temperature field in K. Core dimensions: level, latitude,
        longitude.
    names : dict, optional
        Variable name override.

    Returns
    -------
    xarray.DataArray
        Isentropic density in kg/K/m².
    """
    isob, lat, lon, sg = get_names(names, "isob", "lat", "lon", "sg")
    return xr.apply_ufunc(
        _pvgradient.isentropic_density_isob,
        da_pt[isob],
        da_pt,
        input_core_dims=[
            [isob],
            [isob, lat, lon]
        ],
        output_core_dims=[
            [isob, lat, lon]
        ],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[
            np.float64
        ]
    ).rename(sg)


def interpolate_isob_to_isen(ds, da_isen=None, *, names=None):
    """Interpolation from pressure levels to isentropic levels.

    Either assign the target isentropic levels as a coordinate (isentrope) to
    the input dataset or provide *da_isen*. The pressure coordinate (level)
    must be specified in ascending order, the isentropic levels in descending
    order (both top-down). Out-of-range values are filled with NaN.

    Parameters
    ----------
    ds : xarray.Dataset
        Input data fields. Core dimensions: level, latitude, longitude.
    da_isen : xarray.DataArray, optional
        Override for isentropic levels to interpolate to (in K).
    names : dict, optional
        Variable name override.

    Returns
    -------
    xarray.Dataset
        Dataset with all interpolated fields.
    """
    isen, isob, lat, lon, pt, pres = get_names(names, "isen", "isob", "lat", "lon", "pt", "pres")
    # If isentropic levels are not specified explicitly, take them from the
    # dataset (arity-1 call), otherwise use the supplied levels and their name
    if da_isen is None:
        da_isen = ds[isen]
    else:
        isen = da_isen.name
    # Assemble the varargs list and its in/out configuration for apply_ufunc
    args = [ds[pt]] # potential temperature is always the first in
    nams = [pres] # interpolated pressure is always the first out
    for var in ds.data_vars:
        # Potential temperature is always the first argument, don't include it
        # again in the list of field args
        if var in NAMES["pt"]:
            continue
        # TODO check coords of var (only interpolate those that have the dims)
        args.append(ds[var])
        nams.append(var)
    # Multiple output fields
    outs = xr.apply_ufunc(
        _pvgradient.interpolate_isob_to_isen,
        ds[isob], # isobars (coordinate in)
        da_isen, # isentropes (coordinate out)
        *args, # fields to interpolate
        input_core_dims=[
            [isob], # coordinate in
            [isen], # coordinate out
            *([isob, lat, lon] for _ in args) # fields
        ],
        output_core_dims=[[isen, lat, lon] for _ in nams],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[np.float64 for _ in nams]
    )
    return xr.merge([out.rename(nam) for out, nam in zip(outs, nams)])


def absolute_vorticity(da_u, da_v, *, names=None):
    """Absolute vorticity of a horizontal wind field.

    Coriolis parameter + vertical component of relative vorticity. The planet
    is Earth with angular frequency 7.292e-5 1/s and radius 6371.2 km.

    Parameters
    ----------
    da_u : xarray.DataArray
        Zonal wind component field in m / s. Core dimensions: latitude,
        longitude.
    da_v : xarray.DataArray
        Meridional wind component field in m / s. Core dimensions:
        latitude, longitude.
    names : dict, optional
        Variable name override.

    Returns
    -------
    xarray.DataArray
        Absolute vorticity field in 1 / s.
    """
    lat, lon, av = get_names(names, "lat", "lon", "av")
    return xr.apply_ufunc(
        _pvgradient.absolute_vorticity,
        da_u[lat],
        da_u,
        da_v,
        input_core_dims=[
            [lat],
            [lat, lon],
            [lat, lon]
        ],
        output_core_dims=[
            [lat, lon]
        ],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[np.float64]
    ).rename(av)


def potential_vorticity_isen(da_av, da_sg, *, names=None):
    """Isentropic Ertel PV from potential temperature and isentropic density.

    Potential vorticity q = ζₐ / σ. NaN-values are assigned where
    σ = 0.

    .. note::
        The output is in potential vorticity units. 1 PVU = K m² / s / kg.

    Parameters
    ----------
    da_av : xarray.DataArray
        Absolute vorticity field in 1 / s. Core dimensions: latitude,
        longitude.
    da_sg : xarray.DataArray
        Isentropic density field in kg / K / m². Core dimensions: latitude,
        longitude.
    names : dict, optional
        Variable name override.

    Returns
    -------
    xarray.DataArray
        Potential vorticity in PVU.
    """
    lat, lon, pv = get_names(names, "lat", "lon", "pv")
    return xr.apply_ufunc(
        _pvgradient.potential_vorticity_isen,
        da_av,
        da_sg,
        input_core_dims=[
            [lat, lon],
            [lat, lon]
        ],
        output_core_dims=[
            [lat, lon]
        ],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[np.float64]
    ).rename(pv)


def isob_to_isen_all(ds, da_isen=None, *, vectorize=True, names=None):
    """All-in-one isentropic PV computation from pressure-level inputs.

    Computes pressure (in hPa), isentropic density (in kg / K / m²), absolute
    vorticity (in 1 / s) and Ertel potential vorticity (in PVU) on isentropes
    (in K, descending order) from temperature (in K) and zonal and meridional
    wind components (in m / s) on pressure levels (in hPa, ascending order).

    .. note::
        Using this function is generally faster and easier than going through
        :py:func:`potential_temperature_isob`,
        :py:func:`isentropic_density_isob`,
        :py:func:`interpolate_isob_to_isen`, :py:func:`absolute_vorticity` and
        :py:func:`potential_vorticity_isen` individually.

    Parameters
    ----------
    ds : xarray.Dataset
        Input dataset with temperature and horizontal wind components as well
        as isentropic levels to interpolate to. Core dimensions: level,
        latitude, longitude.
    da_isen : xarray.DataArray, optional
        Override for isentropic levels to interpolate to (in K).
    vectorize : boolean, optional
        Use vectorization of :py:func:`xarray.apply_ufunc`.
    names : dict, optional
        Variable name override.

    Returns
    -------
    xarray.Dataset
        Dataset containing all output variables. Core dimensions: isentrope,
        latitude, longitude.
    """
    isob, isen, lat, lon, t, u, v, pres, sg, av, pv = get_names(
        names, "isob", "isen", "lat", "lon", "t", "u", "v", "pres", "sg", "av", "pv"
    )
    # If isentropic levels are not specified explicitly, take them from the
    # dataset (arity-1 call), otherwise use the supplied levels and their name
    if da_isen is None:
        da_isen = ds[isen]
    else:
        isen = da_isen.name
    # Multiple output fields
    out_names = [pres, u, v, sg, av, pv]
    out = xr.apply_ufunc(
        _pvgradient.isob_to_isen_all,
        ds[isob],
        ds[isen],
        ds[lat],
        ds[t],
        ds[u],
        ds[v],
        #ds[psfc] if psfc in ds else None, # TODO figure out how to input_core_dim this
        input_core_dims=[
            [isob],
            [isen],
            [lat],
            [isob, lat, lon],
            [isob, lat, lon],
            [isob, lat, lon]
        ],
        output_core_dims=[
            [isen, lat, lon],
            [isen, lat, lon],
            [isen, lat, lon],
            [isen, lat, lon],
            [isen, lat, lon],
            [isen, lat, lon]
        ],
        vectorize=vectorize,
        dask="parallelized",
        output_dtypes=[
            np.float64,
            np.float64,
            np.float64,
            np.float64,
            np.float64,
            np.float64
        ]
    )
    return xr.Dataset({ name: field for name, field in zip(out_names, out) })

