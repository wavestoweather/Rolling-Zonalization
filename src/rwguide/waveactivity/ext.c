#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

// Smallest isentropic density value still considered stable and above ground
#define TOL_SG 0.001 // kg / K / m²
#define PVU 1.0e6 // multiply to obtain PV in PVU

// Physical constants
#define RADIUS (6371.2e3) // m

// Trigonometry for degree-valued arguments
#define DEG2RAD(angle) ( (angle) / 180. * M_PI )
#define COSDEG(x) ( cos(DEG2RAD(x)) )
#define SINDEG(x) ( sin(DEG2RAD(x)) )


// Is the gridbox stable based on isentropic density?
inline bool isstable(const double sg) {
    return TOL_SG <= sg;
}

inline bool ispole(const double lat) {
    return fabs(fabs(lat) - 90.) < 0.001;
}

// Spherical domain: line weights
// 2i   values: weight for northern half of arc i
// 2i+1 values: weight for southern half of arc i
inline void boxcount_line_half_weights(
    const double * lat,
    const size_t nlat,
    double * weights
) {
    size_t i = 0;
    // Northern boundary: only southern half of grid cell
    weights[2*i  ] = 0.;
    // Southern half of northernmost grid cell
    weights[2*i+1] = SINDEG(lat[i]) - SINDEG(0.5 * (lat[i] + lat[i+1]));
    // Interior grid cells
    while (++i < nlat - 1) {
        weights[2*i  ] = SINDEG(0.5 * (lat[i-1] + lat[i])) - SINDEG(lat[i]);
        weights[2*i+1] = SINDEG(lat[i]) - SINDEG(0.5 * (lat[i] + lat[i+1]));
    }
    // Northern half of southernmost grid cell
    weights[2*i] = SINDEG(0.5 * (lat[i-1] + lat[i])) - SINDEG(lat[i]);
    // Southern boundary: only northern half of grid cell
    weights[2*i+1] = 0.;
}


void local_wave_activity(
    const double * lat, // [nlat]
    const double * av, // [nlev, nlat, nlon] absolute vorticity in 1/s
    const double * sg, // [nlev, nlat, nlon] isentropic density in kg/K/m²
    const double * pv_bg, // [nlev, nlat, nlon] background-state (zonalized) PV in PVU
    const size_t nlev,
    const size_t nlat,
    const size_t nlon,
    double * lwa // [nlev, nlat, nlon] local wave activity (in m/s)
) {
    const size_t nfld = nlat * nlon;
    // Half-box weights for boxcounting quadrature
    double * weights = malloc(2*nlat * sizeof(double));
    boxcount_line_half_weights(lat, nlat, weights);
    // Precompute full-box weights too
    double * areas = malloc(nlat * sizeof(double));
    for (size_t i = 0; i < nlat; i++) {
        areas[i] = (weights[2*i] + weights[2*i+1]);
    }
    // LWA normalization with a / cos(lat)
    double * normalization = malloc(nlat * sizeof(double));
    for (size_t i = 0; i < nlat; i++) {
        // Avoid div/0 for lat = 90° or -90°
        normalization[i] = ispole(lat[i]) ? 1. : RADIUS / COSDEG(lat[i]);
    }
    // Iterate through levels
    for (size_t h = 0; h < nlev; h++) {
        for (size_t i = 0; i < nlat; i++) {
            for (size_t j = 0; j < nlon; j++) {
                // Contour value
                const double pvc = pv_bg[h*nfld + i*nlon + j] / PVU;
                // Skip if below ground
                if (isnan(pvc)) {
                    lwa[h*nfld + i*nlon + j] = NAN;
                    continue;
                }
                // Initialize LWA accumulator
                lwa[h*nfld + i*nlon + j] = 0.;
                // Local wave activity at the poles is always all zero, so just
                // leave value from initialization
                if (ispole(lat[i])) {
                    continue;
                }
                // Northern integration
                for (size_t k = 0; k < i; k++) {
                    const double avv = av[h*nfld + k*nlon + j];
                    const double sgg = sg[h*nfld + k*nlon + j];
                    if (isstable(sgg) && avv/sgg < pvc) {
                        lwa[h*nfld + i*nlon + j] -= (avv - pvc*sgg) * areas[k];
                    }
                }
                // Current latitude: split gridbox
                const double avv = av[h*nfld + i*nlon + j];
                const double sgg = sg[h*nfld + i*nlon + j];
                if (isstable(sgg)) {
                    if (avv/sgg < pvc) {
                        lwa[h*nfld + i*nlon + j] -= (avv - pvc*sgg) * weights[2*i];
                    } else {
                        lwa[h*nfld + i*nlon + j] += (avv - pvc*sgg) * weights[2*i+1];
                    }
                }
                // Southern integration
                for (size_t k = i+1; k < nlat; k++) {
                    const double avv = av[h*nfld + k*nlon + j];
                    const double sgg = sg[h*nfld + k*nlon + j];
                    if (isstable(sgg) && avv/sgg >= pvc) {
                        lwa[h*nfld + i*nlon + j] += (avv - pvc*sgg) * areas[k];
                    }
                }
                // Normalization with 1/(a*cos(phi))
                lwa[h*nfld + i*nlon + j] *= normalization[i];
            }
        }
    }
    free(weights);
    free(areas);
    free(normalization);
}

