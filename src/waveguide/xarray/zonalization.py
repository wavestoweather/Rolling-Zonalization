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
    return _window(_zonalization.fixed_km_window, da_lon, da_lat, name, **kwargs)

def fixed_deg_window(da_lon, da_lat, name="fixed_deg_window", **kwargs):
    return _window(_zonalization.fixed_deg_window, da_lon, da_lat, name, **kwargs)


def rolling_mean_background(da_area, da_sg, *, vectorize=True, names=None):
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
    lat, lon, pv = _common.get_names(names, "lat", "lon", "pv")
    if da_bg is None:
        da_bg = da_sg
    # If area weights are only a function of latitude, blow-up to 2D field
    if lon not in da_area.coords:
        nlat = da_area.coords[lat].size
        nlon = da_av.coords[lon].size
        area2d = da_area.values.repeat(nlon).reshape((nlat, nlon))
        da_area = xr.DataArray(area2d, coords={
            lat: da_area.coords[lat],
            lon: da_av.coords[lon]
        })
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

