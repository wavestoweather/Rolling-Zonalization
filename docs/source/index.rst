Rossby Waveguide Python Package
===============================

Documentation of Rossby waveguide diagnostics software from `wavestoweather/Rolling-Zonalization <https://github.com/wavestoweather/Rolling-Zonalization>`_.

.. note::
    We are not certain if and how this package will be maintained in the future.
    If you are interested in contributing, don't hesitate to get in touch:

    - https://github.com/wavestoweather/Rolling-Zonalization/issues
    - https://dynmet.ipa.uni-mainz.de/


Contents
--------

.. toctree::
   :maxdepth: 1

   wavenumber
   pvgradient
   zonalization
   visualization
   utils

.. note::
    All submodules have both an xarray and a numpy interface.
    The numpy interface wraps the underlying C-extensions closely.
    The xarray interface wraps the numpy interface, handles coordinates automatically and is `dask-compatible <https://docs.xarray.dev/en/stable/generated/xarray.apply_ufunc.html>`_.
    It is highly recommended to use the xarray-based functions.


Requirements
------------

- C compiler to build extensions
- Python 3 with numpy, scipy, xarray and cffi packages.


Acknowledgements
----------------

This software has been created within the Transregional Collaborative Research Center SFB/TRR 165 "Waves to Weather" funded by the German Science Foundation (DFG). https://www.wavestoweather.de/

