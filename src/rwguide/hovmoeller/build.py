import argparse
import os
import cffi


ffibuilder = cffi.FFI()

code = """
    void mean_along_contour(
        const double xvalue,
        const double * field_x,
        const double * field_y,
        const double * area_weights,
        const double * lat_kernel,
        const double * lon_kernel,
        const size_t nlev,
        const size_t nlat,
        const size_t nlon,
        const size_t nlat_kernel,
        const size_t nlon_kernel,
        double * mean_y
    );
"""
ffibuilder.cdef(code)

ffibuilder.set_source(
    "rwguide.hovmoeller._ext",
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

