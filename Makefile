

.PHONY: docs py-compile py-install


# Python installation

# Version-dependent file ending for compiled C-extensions
PYEXT := $(shell python3-config --extension-suffix)

py-compile: \
		src/waveguide/hovmoeller/_ext$(PYEXT) \
		src/waveguide/pvgradient/_ext$(PYEXT) \
		src/waveguide/wavenumber/_ext$(PYEXT) \
		src/waveguide/zonalization/_ext$(PYEXT)

src/%/_ext$(PYEXT): src/%/build.py src/%/ext.c
	cd src && python3 $*/build.py


# Documentation

# TODO C-extension dependencies
docs: py-compile
	sphinx-build -M html "docs/source" "docs"

