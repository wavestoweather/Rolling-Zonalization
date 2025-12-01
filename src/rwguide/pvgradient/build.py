import argparse
import os
import cffi


ffibuilder = cffi.FFI()

code = """
    int horizontal_gradient(
        const double * lat,
        const double * field,
        const size_t nlev,
        const size_t nlat,
        const size_t nlon,
        double * grad_lon,
        double * grad_lat
    );
    int curl(
        const double * lat,
        const double * u,
        const double * v,
        const size_t nlev,
        const size_t nlat,
        const size_t nlon,
        double * curl
    );
    void norm_grad_log_abs(
        const double * lat,
        const double * pv,
        const double threshold,
        const size_t nlev,
        const size_t nlat,
        const size_t nlon,
        double * out
    );
    int potential_temperature_isob(
        const double * isob,
        const double * t,
        const size_t nisob,
        const size_t nlat,
        const size_t nlon,
        double * th
    );
    int isentropic_density_isob(
        const double * isob,
        const double * th,
        const size_t nisob,
        const size_t nlat,
        const size_t nlon,
        double * sg
    );
    int mask_underground_isob(
        const double * isob,
        const double * psfc,
        double * field,
        const double fill,
        const size_t nisob,
        const size_t nlat,
        const size_t nlon
    );
    int interpolate_isob_to_isen(
        const double * isob,
        const double * isen,
        const double * th,
        const size_t nisob,
        const size_t nisen,
        const size_t nlat,
        const size_t nlon,
        const size_t nvar,
        double * p,
        ...
    );
    int mask_underground_isen(
        const double * p,
        const double * psfc,
        double * field,
        const double fill,
        const size_t nisen,
        const size_t nlat,
        const size_t nlon
    );
    int absolute_vorticity(
        const double * lat,
        const double * u,
        const double * v,
        const size_t nlev,
        const size_t nlat,
        const size_t nlon,
        double * av
    );
    int potential_vorticity_isen(
        const double * av,
        const double * sg,
        const size_t nisen,
        const size_t nlat,
        const size_t nlon,
        double * pv
    );
    void isob_to_isen_all(
        const double * isob,
        const double * isen,
        const double * lat,
        const double * t_isob,
        const double * u_isob,
        const double * v_isob,
        const double * psfc,
        const size_t nvec,
        const size_t nisob,
        const size_t nisen,
        const size_t nlat,
        const size_t nlon,
        double * p_isen,
        double * u_isen,
        double * v_isen,
        double * sg_isen,
        double * av_isen,
        double * pv_isen
    );
"""
ffibuilder.cdef(code)

ffibuilder.set_source(
    "rwguide.pvgradient._ext",
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

