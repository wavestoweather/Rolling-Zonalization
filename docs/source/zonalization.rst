Zonalized Background States
===========================

Zonalization and rolling zonalization of potential vorticity (PV) on isentropes.


Example
-------

>>> from rwguide.xarray import zonalization

Start from a dataset of isentropic density and absolute vorticity on isentropic levels:

>>> isen
<xarray.Dataset>
Dimensions:    (isentrope: 3, latitude: 121, longitude: 240, time: 1460)
Coordinates:
  * isentrope  (isentrope) float64 340.0 330.0 320.0
  * latitude   (latitude) float32 90.0 88.5 87.0 85.5 ... -87.0 -88.5 -90.0
  * longitude  (longitude) float32 -180.0 -178.5 -177.0 ... 175.5 177.0 178.5
  * time       (time) datetime64[ns] 2018-01-01 ... 2018-12-31T18:00:00
Data variables:
    sg         (time, isentrope, latitude, longitude) float64 dask.array<...>
    av         (time, isentrope, latitude, longitude) float64 dask.array<...>

Create a weighting window for rolling zonalization:

>>> win = zonalization.fixed_deg_window(isen["longitude"], isen["latitude"], width=60., window="boxcar")
>>> win
<xarray.DataArray 'fixed_deg_window' (latitude: 121, longitude: 41)>
array(...)
Coordinates:
  * latitude   (latitude) float32 90.0 88.5 87.0 85.5 ... -87.0 -88.5 -90.0
  * longitude  (longitude) float32 0.0 1.5 3.0 4.5 6.0 ... 55.5 57.0 58.5 60.0
Attributes:
    width:    60.0

Compute a rolling-mean isentropic density background:

>>> bgsg = zonalization.rolling_mean_background(win, isen["sg"])
>>> bgsg
<xarray.DataArray 'sg_rm' (time: 1460, isentrope: 3, latitude: 121, longitude: 240)>
dask.array<...>
Coordinates:
  * latitude   (latitude) float32 90.0 88.5 87.0 85.5 ... -87.0 -88.5 -90.0
  * isentrope  (isentrope) float64 340.0 330.0 320.0
  * longitude  (longitude) float32 -180.0 -178.5 -177.0 ... 175.5 177.0 178.5
  * time       (time) datetime64[ns] 2018-01-01 ... 2018-12-31T18:00:00

Run the rolling zonalization:

>>> zonalization.zonalize_rolling(win, isen["av"], isen["sg"], bgsg)
<xarray.DataArray 'pv_rz' (time: 1460, isentrope: 3, latitude: 121, longitude: 240)>
dask.array<...>
Coordinates:
  * latitude   (latitude) float32 90.0 88.5 87.0 85.5 ... -87.0 -88.5 -90.0
  * isentrope  (isentrope) float64 340.0 330.0 320.0
  * longitude  (longitude) float32 -180.0 -178.5 -177.0 ... 175.5 177.0 178.5
  * time       (time) datetime64[ns] 2018-01-01 ... 2018-12-31T18:00:00


Xarray Interface
----------------

.. automodule:: rwguide.xarray.zonalization
    :members:
    :undoc-members:


Numpy Interface
---------------

.. automodule:: rwguide.zonalization
    :members:
    :undoc-members:
    :imported-members:


Weighting Windows
-----------------

The implementation uses a generalized approach to the area weighting of the surface integrals on the sphere evaluated during the zonalization.
This allows for customization of the PV re-arrangement by using arbitrary weighting functions in addition to the normal latitude-dependent cosine-based area weighting.

Regular mass integrals over isentropic density ``σ`` on an isentropic surface of the sphere are given by

    ``∫∫ σ(ϕ,λ) a cos(ϕ) dϕ dλ``,

where ``a cos(ϕ)`` constitutes an area-weighting that accounts for the difference in size of grid boxes of a latitude-longitude grid.
We generalize the weighting to

    ``∫∫ σ(ϕ,λ) w(ϕ,λ) dϕ dλ``,

where ``w`` is now an arbitrary weighting function that can vary both in latitude and longitude.
In the limit of a simple boxcar window, we arrive at basic rolling zonalization where PV is zonalized normally in a rolling window.
But the weighting window can, e.g. be designed with tapered edges and can widen on the lat-lon grid towards the poles to obtain a fixed window width in terms of actual distance.

We provide a set of functions to build customized weighting windows on the sphere based on cosine area weighting and window functions from :py:mod:`scipy.signal`:

.. automodule:: rwguide.zonalization.weighting
    :members:
    :undoc-members:

