# List of years from which to request and process data
YEARS := 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 \
         1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 \
         1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 \
         2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 \
         2019 2020 2021 2022 

# Reanalysis data files, each containing 6-hourly temperature, zonal wind and
# meridional wind on 18 pressure levels at 1.5° horizontal resolution
ERA5_TUV := $(foreach y,$(YEARS),data/ERA5/ERA5-$(y)-tuv-1.5.nc)

# Standard set of isentropes used in the analysis
ISENTROPES := 345 340 335 330 325 320 315#
ISENTROPES_CS := $(subst $(space),$(comma),$(ISENTROPES))# comma-separated



# Version-dependent file ending for compiled C extensions
PYEXT := $(shell python3-config --extension-suffix)

.SECONDARY:
.PHONY: all reanalysis docs py-compile py-install

all: figures/barotropic.pdf


# Rule for creating folders
%/:
	mkdir -p $*


# ERA5 data download
reanalysis: $(ERA5_TUV)

data/ERA5/ERA5-%-tuv-1.5.nc: src/download_ERA5.py | data/ERA5/
	python3 -m src.download_ERA5 $*


# Data processing: basic aggregation
data/mean-isen.nc: src/calculate_means.py $(ERA5_TUV) | data/
	python3 -m src.calculate_means --isen --levels=$(ISENTROPES_CS) $(ERA5_TUV) $@


# Figure 1: barotropic waveguide diagnostics
figures/barotropic.pdf: src/plot_barotropic.py src/common/plotting.py \
		src/waveguide/xarray/wavenumber.py src/waveguide/xarray/pvgradient.py \
		data/mean-isen.nc | figures/
	python3 -m src.plot_barotropic --level=330 --isen --season="winter" data/mean-isen.nc $@


# Python 'waveguide' package
py-install:
	pip install .

py-compile: \
		src/waveguide/hovmoeller/_ext$(PYEXT) \
		src/waveguide/pvgradient/_ext$(PYEXT) \
		src/waveguide/wavenumber/_ext$(PYEXT) \
		src/waveguide/zonalization/_ext$(PYEXT)

src/%/_ext$(PYEXT): src/%/build.py src/%/ext.c
	cd src && python3 $*/build.py

src/waveguide/xarray/%.py: src/waveguide/%/_ext$(PYEXT)
	touch $@


# Documentation
docs: py-compile
	sphinx-build -M html "docs/source" "docs"

