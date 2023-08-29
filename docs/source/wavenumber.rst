Stationary Wavenumber
=====================

Compute the barotropic stationary wavenumber from the zonal wind component.


Example
-------

>>> import xarray as xr
>>> from rwguide.xarray import wavenumber

Loading input data (zonal wind component on pressure levels):

>>> uvt = xr.open_dataset("data/ERA5/ERA5-2018-tuv-1.5.nc", chunks={ "time": 20 })
>>> uvt
<xarray.Dataset>
Dimensions:    (longitude: 240, latitude: 121, level: 18, time: 1460)
Coordinates:
  * longitude  (longitude) float32 -180.0 -178.5 -177.0 ... 175.5 177.0 178.5
  * latitude   (latitude) float32 90.0 88.5 87.0 85.5 ... -87.0 -88.5 -90.0
  * level      (level) int32 50 70 100 150 200 250 ... 600 650 700 750 800 850
  * time       (time) datetime64[ns] 2018-01-01 ... 2018-12-31T18:00:00
Data variables:
    u          (time, level, latitude, longitude) float32 dask.array<...>
    v          (time, level, latitude, longitude) float32 dask.array<...>
    t          (time, level, latitude, longitude) float32 dask.array<...>

Computing stationary wavenumber fields:

>>> wavenumber.stationary_wavenumber(uvt["u"])
<xarray.DataArray 'ks' (time: 1460, level: 18, latitude: 121, longitude: 240)>
dask.array<...>
Coordinates:
  * latitude   (latitude) float32 90.0 88.5 87.0 85.5 ... -87.0 -88.5 -90.0
  * longitude  (longitude) float32 -180.0 -178.5 -177.0 ... 175.5 177.0 178.5
  * level      (level) int32 50 70 100 150 200 250 ... 600 650 700 750 800 850
  * time       (time) datetime64[ns] 2018-01-01 ... 2018-12-31T18:00:00



Xarray Interface
----------------

.. automodule:: rwguide.xarray.wavenumber
    :members:
    :undoc-members:


Numpy Interface
---------------

.. automodule:: rwguide.wavenumber
    :members:
    :undoc-members:
    :imported-members:

