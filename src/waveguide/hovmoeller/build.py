import cffi


if __name__ == "__main__":

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
        "waveguide.hovmoeller._ext",
        code,
        sources=[
            "waveguide/hovmoeller/ext.c",
        ],
        libraries=["m"],
    )
    ffibuilder.compile(verbose=True)

