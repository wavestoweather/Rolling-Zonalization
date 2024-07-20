# Rolling-Zonalization

![Rolling Zonalization example](https://raw.githubusercontent.com/wavestoweather/Rolling-Zonalization/main/docs/source/img/rzon.png)

This repository contains the software to produce all figures of the article: Polster, C., & Wirth, V. (2023). A New Atmospheric Background State to Diagnose Local Waveguidability. *Geophysical Research Letters*, 50, e2023GL106166. https://doi.org/10.1029/2023GL106166.

- Research by [Christopher Polster](https://dynmet.ipa.uni-mainz.de/christopher-polster/) and [Volkmar Wirth](https://dynmet.ipa.uni-mainz.de/volkmar-wirth/).
- Software by Christopher Polster.

[MIT License](LICENSE).

[![DOI](https://zenodo.org/badge/682976838.svg)](https://zenodo.org/badge/latestdoi/682976838)


## The `rwguide` Python Package

An implementation of rolling zonalization is included in this repository.
The software can be installed as a standalone Python package.
Install the package from a clone of the repository:

    $ pip install .

If you are interested in computing (rolling) zonalized background states, start here.
Visit the [package documentation](https://wavestoweather.github.io/Rolling-Zonalization) for more information.
Other content of the repository is primarily included to reproduce the datasets and plots of Polster and Wirth (2023).


## How To Run

Clone this repository:

    $ git clone https://github.com/wavestoweather/Rolling-Zonalization.git


### Software Requirements

- make
- C compiler to build Python extensions
- Python 3 with packages listed in [`requirements.txt`](requirements.txt). The specified versions together with Python 3.9 were used during development, older/newer versions may or may not work. Note that setuptools<70.0.0 is affected by a [security vulnerability](https://github.com/advisories/ghsa-cx63-2mw6-8hw5).


### Data Download, Processing and Plots

    $ make

is configured to download and process all required data and create all figures.
The following steps are carried out:

- Download of ERA5 data (see notes below, `src/download_ERA5.py`).
- Compute climatological means of ERA5 variables and PV on isentropes (`src/calculate_means.py`).
- Plot Figure 1 (`src/plot_barotropic.py`).
- Plot Figure 2 (`src/plot_schematic.py`).
- Compute isentropic Ertel PV from the ERA5 data (`src/calculate_pv.py`) and its zonal wavenumber spectrum of PV (`src/calculate_sprops.py`).
- Compute 14-day rolling mean of PV (`src/calculate_rollmean.py`) and its zonal wavenumber spectrum.
- Compute the zonal wavenumber spectrum of climatological mean PV.
- Compute 90° rolling zonalized PV on 330 K, its climatological mean and zonal wavenumber spectrum and analyze waveguide occurrence.
- Compute 60° rolling zonalized PV on 330 K (`src/calculate_pvrz.py`), its climatological mean and zonal wavenumber spectrum and analyze waveguide occurrence.
- Compute 60° rolling zonalized PV on 345, 340, 335, 325, 320 and 315 K and analyze waveguide occurrence.
- Plot Figure 3 (`src/plot_climatology.py`).
- Plot Figure 4 (`src/plot_episode.py`).

Python extensions are compiled as needed in between.
Use the `--dry-run` option of make to see which commands will be run without executing them.
It is generally not necessary/recommended to parallelize with the `-j` option of make.
The data processing is already parallelized with the default [dask](https://www.dask.org/) scheduler.
Most scripts will provide a progress bar while running.


To only start data downloads without the data processing, use

    $ make reanalysis


### Alternative Input

If you already have a similar dataset containing temperature and horizontal wind components on pressure levels, it should generally be possible to substitute these files.
For reference, the preconfigured ERA5 downloader produces one file per year with content structured as:

    <xarray.Dataset>
    Dimensions:    (longitude: 240, latitude: 121, level: 18, time: 1460)
    Coordinates:
      * longitude  (longitude) float32 -180.0 -178.5 -177.0 ... 175.5 177.0 178.5
      * latitude   (latitude) float32 90.0 88.5 87.0 85.5 ... -87.0 -88.5 -90.0
      * level      (level) int32 50 70 100 150 200 250 ... 600 650 700 750 800 850
      * time       (time) datetime64[ns] 2010-01-01 ... 2010-12-31T18:00:00
    Data variables:
        u          (time, level, latitude, longitude) float32 ...
        v          (time, level, latitude, longitude) float32 ...
        t          (time, level, latitude, longitude) float32 ...

A few changes might be necessary to accomodate different file content structure, grid resolutions, etc. (e.g. by specifying the appropriate file paths in the `Makefile`, adapting the number of timesteps in the window for the rolling temporal mean, ...).


### Output

ERA5 data is placed into `data/ERA5`.
The approximate size of the downloaded dataset is 190 GB.

NetCDF files with intermediate processed fields are written to the `data` directory.
The approximate size of the processed dataset is 35 GB.
File names use prefixes

- `data/PV-*.nc`: potential vorticity,
- `data/PVrm-*.nc`: 14-day rolling-mean potential vorticity,
- `data/PVrz-*.nc`: rolling-zonalized potential vorticity,

and suffixes

- `data/*-mean.nc`: climatological and seasonal means,
- `data/*-sprop.nc`: mean zonal wavenumber spectrum,
- `data/*-occur.nc`: waveguide occurrence frequencies.

Figures are written to the `figures` directory:

- Figure 1: `figures/barotropic.pdf`,
- Figure 2: `figures/schematic.pdf`,
- Figure 3: `figures/climatology.pdf`,
- Figure 4: `figures/episode.pdf`.


### Known Issues

- Scripts occasionally get stuck in the middle of processing while doing seemingly nothing (no significant CPU usage).
  The problem has been difficult to diagnose because the hangups seem to occur at a different stage of the processing every time.
  Terminating a stuck script and re-running `make` can be sufficient to work around the issue.
- Parallel performance can be optimized for a specific machine by setting appropriate chunksizes in the processing scripts, e.g. [here](https://github.com/wavestoweather/Rolling-Zonalization/blob/4e79ec2c6ff63785d2dcff165d172b556c39f7f0/src/calculate_sprops.py#L27)). Additionally, it might be beneficial or necessary to specify `OMP_NUM_THREADS=1` or `OMP_NUM_THREADS=2` to avoid [oversubscription](https://docs.dask.org/en/latest/array-best-practices.html?highlight=OMP_NUM_THREADS#avoid-oversubscribing-threads).


## Acknowledgements

This research project has been carried out within the Transregional Collaborative Research Center SFB/TRR 165 "Waves to Weather" funded by the German Science Foundation (DFG). https://www.wavestoweather.de/

