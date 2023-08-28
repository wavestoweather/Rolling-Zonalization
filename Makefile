# List of years from which to request and process data
YEARS := 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 \
         1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 \
         1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 \
         2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 \
         2019 2020 2021 2022 

# Reanalysis data files, each containing 6-hourly temperature, zonal wind and
# meridional wind on 18 pressure levels at 1.5° horizontal resolution
ERA5_TUV := $(foreach y,$(YEARS),data/ERA5/ERA5-$(y)-tuv-1.5.nc)


# Version-dependent file ending for compiled C extensions
PYEXT := $(shell python3-config --extension-suffix)

.SECONDARY:
.PHONY: all reanalysis docs py-compile py-install

all: figures/barotropic.pdf figures/schematic.pdf figures/climatology.pdf figures/episode.pdf


# Rule for creating folders
%/:
	mkdir -p $*


# Figure 1: barotropic waveguide diagnostics
figures/barotropic.pdf: src/plot_barotropic.py src/common/plotting.py \
		src/waveguide/xarray/wavenumber.py src/waveguide/xarray/pvgradient.py \
		data/mean-isen.nc | figures/
	python3 -m src.plot_barotropic \
			--level=330 \
			--isen \
			--season="winter" \
			data/mean-isen.nc \
			$@

# Figure 2: schematic of rolling zonalization procedure and comparison to 14
#           day mean background state
figures/schematic.pdf: src/plot_schematic.py src/common/plotting.py \
		src/waveguide/xarray/pvgradient.py src/waveguide/xarray/zonalization.py \
		data/ERA5/ERA5-2016-tuv-1.5.nc | figures/
	python3 -m src.plot_schematic \
			--date="2016-12-18T12:00" \
			--scale=60 \
			--taper=0.0 \
			--grad-exclude=0.1 \
			--isentrope=330 \
			--mean-width=14 \
			data/ERA5/ERA5-2016-tuv-1.5.nc \
			$@

# Figure 3: PV spectrum, waveguide occurrence comparison 60° and 90°, vertical
#           waveguide occurrence structure
figures/climatology.pdf: src/plot_climatology.py src/common/plotting.py src/common/seasons.py \
		data/PV-330K-sprop.nc data/PVrm-330K-sprop.nc data/mean-isen-sprop.nc \
		data/PVrz-330K-90deg-occur.nc data/PVrz-330K-90deg-mean.nc data/PVrz-330K-90deg-sprop.nc \
		data/PVrz-330K-60deg-occur.nc data/PVrz-330K-60deg-mean.nc data/PVrz-330K-60deg-sprop.nc \
		data/PVrz-345K-60deg-occur.nc data/PVrz-340K-60deg-occur.nc data/PVrz-335K-60deg-occur.nc \
		data/PVrz-325K-60deg-occur.nc data/PVrz-320K-60deg-occur.nc data/PVrz-315K-60deg-occur.nc \
		| figures/
	python3 -m src.plot_climatology \
			--scale=60 \
			--levels="345,340,335,330,325,320,315" \
			--level-cmp=330 \
			--scale-cmp=90 \
			--season=winter \
			--threshold="1.2e-6" \
			$@

# Figure 4: Dec 2016 and Jan 2018 wave propagation episodes
figures/episode.pdf: src/plot_episode.py src/common/plotting.py \
		src/waveguide/xarray/pvgradient.py src/waveguide/xarray/hovmoeller.py \
		data/ERA5/ERA5-2016-tuv-1.5.nc data/ERA5/ERA5-2018-tuv-1.5.nc \
		data/PVrz-330K-60deg.nc | figures/
	python3 -m src.plot_episode \
			--scale=60 \
			--grad-exclude=0.1 \
			--isentrope=330 \
			$@


# ERA5 data download
reanalysis: $(ERA5_TUV)

data/ERA5/ERA5-%-tuv-1.5.nc: src/download_ERA5.py | data/ERA5/
	python3 -m src.download_ERA5 $*


# Data processing: basic aggregation
data/mean-isen.nc: src/calculate_means.py $(ERA5_TUV) | data/
	python3 -m src.calculate_means --isen --levels="345,340,335,330,325,320,315" $(ERA5_TUV) $@

data/%-mean.nc: data/%.nc src/calculate_means.py
	python3 -m src.calculate_means $< $@

# Data processing: spatial properties
data/%-sprop.nc: src/calculate_sprops.py data/%.nc
	python3 -m src.calculate_sprops data/$*.nc $@

# Data processing: waveguide occurrence frequency
data/%-occur.nc: data/%.nc src/calculate_occurrence.py
	python3 -m src.calculate_occurrence --threshold="0.8e-6,1.0e-6,1.2e-6,1.4e-6" --grad-exclude=0.1 $< $@

# Background state: just PV
data/PV-%K.nc: src/calculate_pv.py $(ERA5_TUV)
	python3 -m src.calculate_pv --isentropes=$* $(ERA5_TUV) $@

# Background state: rolling-temporal-mean PV (57 = (14 * 24h / 6h) + 1 = 2 weeks of 6-hourly data)
data/PVrm-%K.nc: src/calculate_rollmean.py data/PV-%K.nc
	python3 -m src.calculate_rollmean --length=57 --window=boxcar data/PV-$*K.nc $@

# Background state: rolling zonalized PV (standard window)
data/PVrz-%K-60deg.nc: src/calculate_pvrz.py $(ERA5_TUV)
	python3 -m src.calculate_pvrz --scale=60 --taper=0.0 --isentropes=$* $(ERA5_TUV) $@

# Background state: additional rule for different window widths on 330K
data/PVrz-330K-%deg.nc: src/calculate_pvrz.py $(ERA5_TUV)
	python3 -m src.calculate_pvrz --scale=$* --taper=0.0 --isentropes=330 $(ERA5_TUV) $@


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

