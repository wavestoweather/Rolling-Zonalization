#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Physical constants
#define RADIUS (6371.2e3) // m
#define OMEGA  (7.292e-5) // 1/s

// Trigonometry for degree-valued arguments
#define DEG2RAD(angle) ( (angle) / 180. * M_PI )
#define COSDEG(x) ( cos(DEG2RAD(x)) )

// Variables in the functions defined here follow a common naming convention.
// Define a few macros as shorthands for common computations.
// 3D index macros:
#define IDXN ( h*nfield + (i-1)*nlon + j )
#define IDX0 ( h*nfield +  i   *nlon + j )
#define IDXS ( h*nfield + (i+1)*nlon + j )
// Cosine macros for various (half-)latitudes:
#define COSN  ( COSDEG(lat[i-1]) )
#define COSNH ( COSDEG(0.5 * (lat[i] + lat[i-1])) )
#define COS0  ( COSDEG(lat[i]) )
#define COSSH ( COSDEG(0.5 * (lat[i] + lat[i+1])) )
#define COSS  ( COSDEG(lat[i+1]) )


int stationary_wavenumber_squared(
    double * lat,
    double * u,
    size_t nlev,
    size_t nlat,
    size_t nlon,
    double * ks2
) {
    size_t nfield = nlat * nlon;
    double dphi = DEG2RAD(lat[1] - lat[0]);
    // Pre-compute coefficients for finite-difference formula. Only needed in
    // the interior domain since Ks² is set 0 at the poles.
    double * a = malloc(nlat * sizeof(nlat));
    double * b = malloc(nlat * sizeof(nlat));
    double * c = malloc(nlat * sizeof(nlat));
    double * d = malloc(nlat * sizeof(nlat));
    for (size_t i = 1; i < nlat-1; i++) {
        a[i] = 2. * OMEGA * RADIUS * COS0 * COS0 * COS0;
        b[i] = - COS0 * COS0 * COSN / dphi / dphi / COSNH;
        c[i] =   COS0 * COS0 * COS0 / dphi / dphi * (1./COSNH + 1./COSSH);
        d[i] = - COS0 * COS0 * COSS / dphi / dphi / COSSH;
    }
    for (size_t h = 0; h < nlev; h++) {
        for (size_t j = 0; j < nlon; ++j) {
            ks2[h*nfield + j] = 0.0; // North Pole
        }
        for (size_t i = 1; i < nlat-1; i++) {
            for (size_t j = 0; j < nlon; j++) {
                ks2[IDX0] = (a[i] + b[i]*u[IDXN] + c[i]*u[IDX0] + d[i]*u[IDXS]) / u[IDX0];
            }
        }
        for (size_t j = 0; j < nlon; ++j) {
            ks2[h*nfield + (nlat-1)*nlon + j] = 0.0; // South Pole
        }
    }
    free(a);
    free(b);
    free(c);
    free(d);
    return 0;
}


int stationary_wavenumber(
    double * lat,
    double * u,
    size_t nlev,
    size_t nlat,
    size_t nlon,
    double * ks
) {
    int status = stationary_wavenumber_squared(lat, u, nlev, nlat, nlon, ks);
    if (status != 0) {
        return status;
    }
    for (size_t i = 0; i < nlev * nlat * nlon; ++i) {
        ks[i] = (ks[i] >= 0) ? sqrt(ks[i]) : -sqrt(-ks[i]);
    }
    return 0;
}

