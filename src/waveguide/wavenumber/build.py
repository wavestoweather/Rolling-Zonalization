import cffi


if __name__ == "__main__":

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
            "waveguide/wavenumber/ext.c",
        ],
        libraries=["m"],
    )
    ffibuilder.compile(verbose=True)

