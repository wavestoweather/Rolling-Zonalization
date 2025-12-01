import argparse
import os
import cffi


ffibuilder = cffi.FFI()

code = """
    int stationary_wavenumber_squared(double * lat, double * u, size_t nlev, size_t nlat, size_t nlon, double * ks2);
    int stationary_wavenumber(double * lat, double * u, size_t nlev, size_t nlat, size_t nlon, double * ks);
"""
ffibuilder.cdef(code)

ffibuilder.set_source(
    "rwguide.wavenumber._ext",
    code,
    sources=[
        os.path.join(os.path.dirname(os.path.relpath(__file__)), "ext.c")
    ],
    libraries=["m"],
)

parser = argparse.ArgumentParser()
parser.add_argument("--emit-c-code", type=str, default=None)

if __name__ == "__main__":
    args = parser.parse_args()

    if args.emit_c_code is not None:
        ffibuilder.emit_c_code(args.emit_c_code)
    else:
        ffibuilder.compile(verbose=True)
