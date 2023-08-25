import os
import cffi


ffibuilder = cffi.FFI()

code = """
    int stationary_wavenumber_squared(double * lat, double * u, size_t nlev, size_t nlat, size_t nlon, double * ks2);
    int stationary_wavenumber(double * lat, double * u, size_t nlev, size_t nlat, size_t nlon, double * ks);
"""
ffibuilder.cdef(code)

ffibuilder.set_source(
    "waveguide.wavenumber._ext",
    code,
    sources=[
        os.path.join(os.path.dirname(os.path.relpath(__file__)), "ext.c")
    ],
    libraries=["m"],
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)

