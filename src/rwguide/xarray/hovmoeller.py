import numpy as np
import xarray as xr

from . import common as _common
from .. import hovmoeller as _hovmoeller


def mean_along_contour(xvalue, da_x, da_y, da_area=None, lat_kernel=None, lon_kernel=None, *, vectorize=True, names=None):
    """Kernel-weighted mean along a contour.

    .. note::
        This function can be used to create a refined Hovmöller diagram as
        described by `Martius et al. (2006)`_ along a contour of a rolling
        zonalized state, since almost all PV contours in a rolling zonalized
        state are circumglobal and occur exactly once. It should not be used to
        create refined Hovmöller diagrams of arbitrary input fields, where
        contours do not fulfil these properties. No checks on the followed
        contour are performed.

    .. _Martius et al. (2006): https://dx.doi.org/10.1111/j.1600-0870.2006.00172.x

    Parameters
    ----------
    xvalue : number
        Value of the contour to follow.
    da_x : xarray.DataArray
        Field containing the contour to follow. Core dimensions: latitude,
        longitude.
    da_y : xarray.DataArray
        Field of which the mean is taken. Core dimensions: latitude, longitude.
    da_area : xarray.DataArray, optional
        Area weights for a weighted mean on the sphere. Core dimension:
        latitude. Weights are normalized during the computation. If *None*,
        ``cos(ϕ)``-weights are used.
    lat_kernel : numpy.ndarray, optional
        Latitude weighting kernel about the contour. Size in gridpoints.
        Weights are normalized during the computation.
    lon_kernel : numpy.ndarray, optional
        Longitude weighting kernel about the contour. Size in gridpoints.
        Weights are normalized during the computation.
    vectorize : boolean, optional
        Use vectorization of :py:func:`xarray.apply_ufunc`.
    names : dict, optional
        Variable name override.

    Returns
    -------
    xarray.DataArray
        Mean of *da_y* along contour *xvalue* in *da_x*.

    """
    lat, lon = _common.get_names(names, "lat", "lon")
    if da_area is None:
        da_area = np.cos(np.deg2rad(da_x.coords[lat]))
    return xr.apply_ufunc(
        _hovmoeller.mean_along_contour,
        xvalue,
        da_x,
        da_y,
        da_area,
        kwargs={
            "lat_kernel": lat_kernel,
            "lon_kernel": lon_kernel
        },
        input_core_dims=[
            [],
            [lat, lon],
            [lat, lon],
            [lat]
        ],
        output_core_dims=[
            [lon]
        ],
        vectorize=vectorize,
        dask="parallelized",
        output_dtypes=[
            np.float64
        ]
    ).rename(da_y.name)

