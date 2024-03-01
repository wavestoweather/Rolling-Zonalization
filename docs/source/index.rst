Rossby Waveguide Python Package
===============================

Documentation of Rossby waveguide diagnostics software from `wavestoweather/Rolling-Zonalization <https://github.com/wavestoweather/Rolling-Zonalization>`_.

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.8300732.svg
  :target: https://doi.org/10.5281/zenodo.8300732


Contents
--------

.. toctree::
   :maxdepth: 1

   wavenumber
   pvgradient
   zonalization
   waveactivity
   visualization
   utils

.. note::
    All submodules have both an xarray and a numpy interface.
    The numpy interface wraps the underlying C-extensions closely.
    The xarray interface wraps the numpy interface, handles coordinates automatically and is `dask-compatible <https://docs.xarray.dev/en/stable/generated/xarray.apply_ufunc.html>`_.
    It is highly recommended to use the xarray-based functions.

    Many functions offer to use the vectorization of :py:func:`xarray.apply_ufunc`.
    This vectorization mode can be slower than the vectorization built into the C-extensions but often uses much less memory.


Setup
-----

**Requirements:**

- C compiler to build extensions
- Python 3 with numpy, scipy, xarray, cffi and setuptools packages.

**Local install:**

Clone or download the `Rolling-Zonalization <https://github.com/wavestoweather/Rolling-Zonalization>`_ repository and run

.. code-block:: bash

    $ pip install .

from the root of the repository.

News
----

**Release:** Version 1.2 (27 Dec 2023)

- New :py:mod:`rwguide.waveactivity` submodule with a basic implementation of finite-amplitude local wave activity.
- Added `da_isen` argument for :py:func:`rwguide.xarray.pvgradient.isob_to_isen_all`.
- Improved documentation.


**Release:** Version 1.1 (17 Nov 2023)

- Fixed an issue about repeated calls to ``free`` in :py:func:`rwguide.pvgradient.horizontal_gradient`.
- Improved documentation.


Acknowledgements
----------------

This software has been created within the Transregional Collaborative Research Center SFB/TRR 165 "Waves to Weather" funded by the German Science Foundation (DFG). https://www.wavestoweather.de/

