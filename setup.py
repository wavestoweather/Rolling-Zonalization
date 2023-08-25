from setuptools import setup, find_packages

setup(
    name="waveguide",
    description="Waveguide analysis",
    version="1.0.0",
    author="Christopher Polster",
    url="https://github.com/wavestoweather/Rolling-Zonalization",
    packages=find_packages(where="src"),
    install_requires=[
        "cffi",
        "numpy",
        "scipy",
        "xarray"
    ],
    package_dir={
        "": "src",
    },
    cffi_modules=[
        "src/waveguide/hovmoeller/build.py:ffibuilder",
        "src/waveguide/pvgradient/build.py:ffibuilder",
        "src/waveguide/wavenumber/build.py:ffibuilder",
        "src/waveguide/zonalization/build.py:ffibuilder"
    ]
)

