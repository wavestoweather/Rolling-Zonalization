Isentropic PV Gradient
======================

Interpolate to isentropic levels and compute potential vorticity (PV) and its gradient as a waveguide diagnostic.


Example
-------

>>> import xarray as xr
>>> from rwguide.xarray import pvgradient

Loading input data (temperature, wind components on pressure levels):

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

Interpolating to isentropes and computing PV and associated variables:

>>> isen = pvgradient.isob_to_isen_all(uvt.assign_coords({ "isentrope": [340., 330., 320.]}))
>>> isen
<xarray.Dataset>
Dimensions:    (isentrope: 3, latitude: 121, longitude: 240, time: 1460)
Coordinates:
  * isentrope  (isentrope) float64 340.0 330.0 320.0
  * latitude   (latitude) float32 90.0 88.5 87.0 85.5 ... -87.0 -88.5 -90.0
  * longitude  (longitude) float32 -180.0 -178.5 -177.0 ... 175.5 177.0 178.5
  * time       (time) datetime64[ns] 2018-01-01 ... 2018-12-31T18:00:00
Data variables:
    pres       (time, isentrope, latitude, longitude) float64 dask.array<...>
    u          (time, isentrope, latitude, longitude) float64 dask.array<...>
    v          (time, isentrope, latitude, longitude) float64 dask.array<...>
    sg         (time, isentrope, latitude, longitude) float64 dask.array<...>
    av         (time, isentrope, latitude, longitude) float64 dask.array<...>
    pv         (time, isentrope, latitude, longitude) float64 dask.array<...>

Computing the magnitude of the gradient of the logarithm of PV as
a waveguidability proxy:

>>> pvgradient.norm_grad_log_abs(isen["pv"])
<xarray.DataArray 'ngl_pv' (time: 1460, isentrope: 3, latitude: 121, longitude: 240)>
dask.array<...>
Coordinates:
  * latitude   (latitude) float32 90.0 88.5 87.0 85.5 ... -87.0 -88.5 -90.0
  * isentrope  (isentrope) float64 340.0 330.0 320.0
  * longitude  (longitude) float32 -180.0 -178.5 -177.0 ... 175.5 177.0 178.5
  * time       (time) datetime64[ns] 2018-01-01 ... 2018-12-31T18:00:00


Xarray Interface
----------------

.. automodule:: rwguide.xarray.pvgradient
    :members:
    :undoc-members:


Numpy Interface
---------------

.. automodule:: rwguide.pvgradient
    :members:
    :undoc-members:
    :imported-members:

