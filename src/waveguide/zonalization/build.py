import cffi


if __name__ == "__main__":

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
        "waveguide.zonalization._ext",
        code,
        sources=[
            "waveguide/zonalization/ext.c",
        ],
        libraries=["m"],
    )
    ffibuilder.compile(verbose=True)

