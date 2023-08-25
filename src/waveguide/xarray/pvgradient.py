import numpy as np
import xarray as xr

from .common import NAMES, get_names
from .. import pvgradient as _pvgradient


def horizontal_gradient(da, *, vectorize=True, names=None):
    """Horizontal gradient ``∇·`` on the sphere.

    The sphere has the dimensions of Earth (``6371.2 km`` radius).

    Parameters
    ----------
    da : xarray.DataArray
        Input field.
    vectorize : boolean, optional
        Use xarray's vectorization, which can be slower than the built-in
        vectorization of the C-extension but often uses less memory.
    names : dict, optional
        Override for automatic coordinate name recognition.

    Returns
    -------
    Tuple[xarray.DataArray, xarray.DataArray]
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
    """Vertical component of rotation ``∇ ⨯ ·`` on the sphere.

    The sphere has the dimensions of Earth (``6371.2 km`` radius).

    Parameters
    ----------
    da_u : xarray.DataArray
        Input zonal vector component.
    da_v : xarray.DataArray
        Input meridional vector component.
    vectorize : boolean, optional
        Use xarray's vectorization, which can be slower than the built-in
        vectorization of the C-extension but often uses less memory.
    names : dict, optional
        Override for automatic coordinate name recognition.

    Returns
    -------
    xarray.DataArray
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
    """Compute ``|∇(log(|·|))|₂``, as a waveguide diagnostic of isentropic PV.

    To first, order ``∇(log(|IPV|)) ≈ 1/f₀ ∇QGPV ~ ∇²U`` `(Martius et al.
    2009)`_. 

    .. _(Martius et al. 2009): https://dx.doi.org/10.1175/2009JAS2995.1

    Parameters
    ----------
    da : xarray.DataArray
        Input field. Must have latitude and longitude dimensions.
    threshold : number, optional
        Threshold underneath which values of ``|·|`` are ignored, to deal with
        regions where the input field approaches zero and the logarithm tends
        toward negative infinity.
    vectorize : boolean, optional
        Use xarray's vectorization, which can be slower than the built-in
        vectorization of the C-extension but often uses less memory.
    names : dict, optional
        Override for automatic coordinate name recognition.

    Returns
    -------
    xarray.DataArray
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


def isentropic_density_isob(da_th, *, names=None):
    isob, lat, lon, sg = get_names(names, "isob", "lat", "lon", "sg")
    return xr.apply_ufunc(
        _pvgradient.isentropic_density_isob,
        da_th[isob],
        da_th,
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


def isob_to_isen_all(ds, *, vectorize=True, names=None):
    """Isentropic PV from isobaric horizontal wind and temperature.

    Compute pressure, isentropic density, absolute vorticity and Ertel
    potential vorticity on isentropes.

    Parameters
    ----------
    ds : xarray.Dataset
        Input dataset. Must contain zonal wind (m/s), meridional wind (m/s) and
        temperature (K) on pressure levels (hPa) as well as isentropic levels
        (K) to interpolate to.
    vectorize : boolean, optional
        Use xarray's vectorization, which can be slower than the built-in
        vectorization of the C-extension but often uses less memory.
    names : dict, optional
        Override for automatic coordinate name recognition.

    Returns
    -------
    xarray.Dataset
    """
    isob, isen, lat, lon, t, u, v, pres, sg, av, pv = get_names(
        names, "isob", "isen", "lat", "lon", "t", "u", "v", "pres", "sg", "av", "pv"
    )
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

