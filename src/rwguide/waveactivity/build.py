import argparse
import os
import cffi


ffibuilder = cffi.FFI()

code = """
    void local_wave_activity(
        const double * lat,
        const double * av,
        const double * sg,
        const double * pv_bg,
        const size_t nlev,
        const size_t nlat,
        const size_t nlon,
        double * lwa
    );
"""
ffibuilder.cdef(code)

ffibuilder.set_source(
    "rwguide.waveactivity._ext",
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
