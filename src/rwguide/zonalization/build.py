import argparse
import os
import cffi


ffibuilder = cffi.FFI()

code = """
    int rolling_mean_background(
        const double * window,
        const double * field,
        const size_t nlev,
        const size_t nlat,
        const size_t nlon,
        const size_t nwin,
        double * out
    );
    void zonalize(
        const double * area,
        const double * av,
        const double * sg,
        const double * bg,
        const size_t nh_last,
        const size_t sh_first,
        const size_t nlev,
        const size_t nlat,
        const size_t nlon,
        double * pv_zonalized
    );
    void zonalize_rolling(
        const double * area,
        const double * av,
        const double * sg,
        const double * bg,
        const size_t nh_last,
        const size_t sh_first,
        const size_t nlev,
        const size_t nlat,
        const size_t nlon,
        const size_t nwin,
        double * pv_zonalized
    );
"""
ffibuilder.cdef(code)

ffibuilder.set_source(
    "rwguide.zonalization._ext",
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
