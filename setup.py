from setuptools import setup, find_packages

setup(
    name="rwguide",
    description="Rossby waveguide analysis",
    version="1.2.1",
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
        "src/rwguide/hovmoeller/build.py:ffibuilder",
        "src/rwguide/pvgradient/build.py:ffibuilder",
        "src/rwguide/wavenumber/build.py:ffibuilder",
        "src/rwguide/zonalization/build.py:ffibuilder"
    ]
)

