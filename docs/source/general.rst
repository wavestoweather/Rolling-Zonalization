General Information
===================


Interfaces
----------

All submodules have both an xarray and a numpy interface.
The numpy interface wraps the underlying C-extensions closely.
The xarray interface wraps the numpy interface, handles coordinates automatically and is `dask-compatible <https://docs.xarray.dev/en/stable/generated/xarray.apply_ufunc.html>`_.
It is highly recommended to use the xarray-based functions.


Coordinates and Units
---------------------

- Input data must be provided on a **regular longitude-latitude grid**.
  Grid spacing is generally allowed to differ between the longitude and latitude dimensions but must be equally spaced within each dimension.
  E.g., a grid with resolution 1°×2° should work, but a Gaussian grid will not.
- **Latitude coordinates** must start in the north.
  The North Pole is located at +90°, the South Pole at -90°.
- **Vertical coordinates** must be specified top-down if the vertical coordinate is a core dimension.
  For a pressure coordinate that means values must be increasing (ascending order), for an isentropic coordinate values must be decreasing (descending order).

The following names for data coordinates are recognized by the xarray interface:

.. list-table::
    :align: left
    :width: 100%
    :widths: 10 35 10 35 10
    :header-rows: 1

    * - Name
      - Description
      - Unit
      - Notes
      - Override

    * - ``latitude``
      - latitude
      - °
      - N to S
      - ``lat``

    * - ``longitude``
      - longitude
      - °
      - 
      - ``lon``

    * - ``level``
      - pressure (isobaric levels)
      - hPa
      - ascending if core dim
      - ``isob``

    * - ``isentrope``
      - potential temperature (isentropic levels)
      - K
      - descending if core dim
      - ``isen``


The following names for data fields are recognized by the xarray interface:


.. list-table::
    :align: left
    :width: 100%
    :widths: 10 35 10 35 10
    :header-rows: 1
 
    * - Name
      - Description
      - Unit
      - Notes
      - Override
 
    * - ``u``
      - zonal wind
      - m/s
      - positive towards east
      - ``u``
 
    * - ``v``
      - meridional wind
      - m/s
      - positive towards north
      - ``v``
 
    * - ``t``
      - temperature
      - K
      - 
      - ``t``
 
    * - ``pres``
      - pressure
      - hPa
      - 
      - ``pres``
 
    * - ``pt``
      - potential temperature
      - K
      - 
      - ``pt``
 
    * - ``sg``
      - isentropic density
      - kg/K/m²
      - 
      - ``sg``
 
    * - ``vo``
      - relative vorticity
      - 1/s
      - 
      - ``vo``
 
    * - ``av``
      - absolute vorticity
      - 1/s
      - 
      - ``av``
 
    * - ``pv``
      - Ertel potential vorticity
      - PVU
      - 
      - ``pv``
 
    * - ``ks``
      - stationary wavenumber
      - unitless
      - 
      - ``ks``

Overrides can be provided for nonconforming data in any function of the xarray interface that accepts a `names` argument.
E.g., if a pressure coordinate is called ´´isobarticInhPa´´ instead of ´´level´´, use an override (based on the override name from the tables above):

>>> rwguide.xarray.pvgradient.potential_temperature_isob(..., names={ "iosb": "isobaricInhPa" })


Vectorization
-------------

Many functions offer to use the vectorization of :py:func:`xarray.apply_ufunc`.
This vectorization mode can be slower than the vectorization built into the C-extensions but often uses much less memory.

