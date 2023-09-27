#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdbool.h>

// Physical constants
#define P_REF (1000.) // hPa
#define KAPPA (0.286) // dimless
#define GRAVITY (9.80655) // m/s²
#define RADIUS (6371200.) // m
#define OMEGA (7.292e-5) // 1/s
#define PVU 1.0e6 // multiply to obtain PV in PVU
// Smallest isentropic density value still considered stable and above ground
#define TOL_SG 0.001 // kg / K / m²

// Assign only if the value is not NaN
#define ASSIGN_PRESERVE_NAN(var, val) if (!isnan(var)) { var = val; }


// Degree to radians
inline double deg2rad(const double angle) {
    return angle / 180. * M_PI;
}

// Distance between equidistant gridpoints on a closed interval
inline double equidistant_delta(const double min, const double max, const size_t n) {
    return (max - min) / ((double) (n - 1));
}

// Is x contained in the interval [x0, x1) if x0 < x1 or [x1, x0) if x1 < x0?
inline bool isin(const double x, const double x0, const double x1) {
    return (x0 <= x) ? (x < x1) : (x1 <= x);
}

// Is the latitude (in degrees) a pole (90°N or 90°S)?
inline bool ispole(const double latitude) {
    return fabs(fabs(latitude) - 90.) < 0.001;
}

// Is the gridbox stable based on isentropic density?
inline bool isstable(const double sg) {
    return TOL_SG <= sg;
}

// Is the gridbox unstable based on isentropic density?
inline bool isunstable(const double sg) {
    return sg < TOL_SG;
}


// Mean value, skipping NaNs. NaN is returned only if all values are NaN.
inline double nanmean(const double * array, const size_t n) {
    size_t nan = 0;
    double sum = 0.;
    for (size_t i = 0; i < n; i++) {
        if (isnan(array[i])) {
            nan += 1;
        } else {
            sum += array[i];
        }
    }
    return (n == nan) ? NAN : (sum / ((double) (n - nan)));
}

// First Derivative finite-difference kernel for equidistant grid. Uses
// a second-order three-point stencil if there are three valid gridpoints,
// otherwise a first-order two-point stencil.
inline double fd_kernel_1_eq32(
    const double y0,
    const double y1,
    const double y2,
    const double delta
) {
    if (isnan(y1)) {
        return NAN;
    } else if (isnan(y0)) {
        return (y2 - y1) / delta;
    } else if (isnan(y2)) {
        return (y1 - y0) / delta;
    } else {
        return (y2 - y0) / (2. * delta);
    }
}

// Horizontal gradient of a scalar field
int horizontal_gradient(
    const double * lat, // [nlat] in: latitude grid
    const double * field, // [nlev, nlat, nlon] in: scalar field
    const size_t nlev,
    const size_t nlat,
    const size_t nlon,
    double * grad_lon, // [nlev, nlat, nlon] out: zonal gradient
    double * grad_lat  // [nlev, nlat, nlon] out: meridional gradient
) {
    const size_t nfield = nlat * nlon;
    // Grid spacings: assume equidistant grid
    // Latitudes can be pole to pole or only (part of) a hemisphere
    const double dphi = deg2rad(equidistant_delta(lat[0], lat[nlat-1], nlat));
    // Longitudes always have to contain the entire 360°, periodic
    const double dlambda = deg2rad(equidistant_delta(0., 360., nlon+1));
    // Precompute coefficients
    const double latcoef = 1. / RADIUS;
    double * loncoef = malloc(nlat * sizeof(double));
    for (size_t i = 0; i < nlat; i++) {
        loncoef[i] = 1. / RADIUS / cos(deg2rad(lat[i]));
    }
    // Vertical levels are independent
    for (size_t h = 0; h < nlev; h++) {
        // Longitudinal direction: d/dlambda term
        for (size_t i = 0; i < nlat; i++) {
            // Western boundary
            grad_lon[h*nfield + i*nlon] = loncoef[i] * fd_kernel_1_eq32(
                field[h*nfield + i*nlon + nlon-1], // wrapped around periodic boundary
                field[h*nfield + i*nlon         ],
                field[h*nfield + i*nlon + 1     ],
                dlambda
            );
            // Interior domain
            for (size_t j = 1; j < nlon-1; j++) {
                grad_lon[h*nfield + i*nlon + j] = loncoef[i] * fd_kernel_1_eq32(
                    field[h*nfield + i*nlon + j-1],
                    field[h*nfield + i*nlon + j  ],
                    field[h*nfield + i*nlon + j+1],
                    dlambda
                );
            }
            // Eastern boundary
            grad_lon[h*nfield + i*nlon + nlon-1] = loncoef[i] * fd_kernel_1_eq32(
                field[h*nfield + i*nlon + nlon-2],
                field[h*nfield + i*nlon + nlon-1],
                field[h*nfield + i*nlon         ], // wrapped around periodic boundary
                dlambda
            );
        }

        // Meridional direction: d/dphi term
        // Interior domain
        for (size_t i = 1; i < nlat-1; i++) {
            for (size_t j = 0; j < nlon; j++) {
                grad_lat[h*nfield + i*nlon + j] = latcoef * fd_kernel_1_eq32(
                    field[h*nfield + (i-1)*nlon + j],
                    field[h*nfield + (i  )*nlon + j],
                    field[h*nfield + (i+1)*nlon + j],
                    dphi
                );
            }
        }
        // Northern boundary if it is the North Pole: set all values to 0.
        if (ispole(lat[0])) {
            for (size_t j = 0; j < nlon; j++) {
                ASSIGN_PRESERVE_NAN(grad_lat[h*nfield + j], 0.);
                ASSIGN_PRESERVE_NAN(grad_lon[h*nfield + j], 0.);
            }
        // Northern boundary otherwise: use one-sided difference
        } else {
            for (size_t j = 0; j < nlon; j++) {
                grad_lat[h*nfield + j] = latcoef * fd_kernel_1_eq32(
                    NAN, // force one-sided difference
                    field[h*nfield        + j],
                    field[h*nfield + nlon + j],
                    dphi
                );
            }
        }
        // Southern boundary if it is the South Pole: set all values to 0.
        if (ispole(lat[nlat-1])) {
            for (size_t j = 0; j < nlon; j++) {
                ASSIGN_PRESERVE_NAN(grad_lat[h*nfield + (nlat-1)*nlon + j], 0.);
                ASSIGN_PRESERVE_NAN(grad_lon[h*nfield + (nlat-1)*nlon + j], 0.);
            }
        // Southern boundary otherwise: use one-sided difference
        } else {
            for (size_t j = 0; j < nlon; j++) {
                grad_lat[h*nfield + (nlat-1)*nlon + j] = latcoef * fd_kernel_1_eq32(
                    field[h*nfield + (nlat-2)*nlon + j],
                    field[h*nfield + (nlat-1)*nlon + j],
                    NAN, // force one-sided difference
                    dphi
                );
            }
        }
    }
    free(loncoef);
    return 0;
}

// Vertical component of the rotation of a vector field
int curl(
    const double * lat, // [lat] in: latitude grid
    const double * u, // [nlev, nlat, nlon] in: zonal component
    const double * v, // [nlev, nlat, nlon] in: meridional component
    const size_t nlev,
    const size_t nlat,
    const size_t nlon,
    double * curl // [nlev, nlat, nlon] out: curl
) {
    const size_t nfield = nlat * nlon;
    // Grid spacings: assume equidistant grid
    // Latitudes can be pole to pole or only (part of) a hemisphere
    const double dphi = (lat[nlat-1] - lat[0]) / ((double) (nlat - 1)) / 180. * M_PI;
    // Longitudes always have to contain the entire 360°
    const double dlambda = 2. * M_PI / ((double) nlon);
    // Precompute coefficients
    double * coslat = malloc(nlat * sizeof(double)); // cosine of latitude
    double * curlcoef = malloc(nlat * sizeof(double)); // factor for curl in spherical coordinates
    for (size_t i = 0; i < nlat; ++i) {
        coslat[i] = cos(deg2rad(lat[i]));
        curlcoef[i] = 1. / RADIUS / cos(deg2rad(lat[i]));
    }
    // Vertical levels are independent
    for (size_t h = 0; h < nlev; h++) {
        // Relative vorticity, dv/dlambda term, periodic boundaries
        for (size_t i = 0; i < nlat; i++) {
            // Western boundary
            curl[h*nfield + i*nlon] = curlcoef[i] * fd_kernel_1_eq32(
                v[h*nfield + i*nlon + nlon-1], // wrapped around periodic boundary
                v[h*nfield + i*nlon         ],
                v[h*nfield + i*nlon + 1     ],
                dlambda
            );
            // Interior domain
            for (size_t j = 1; j < nlon-1; j++) {
                curl[h*nfield + i*nlon + j] = curlcoef[i] * fd_kernel_1_eq32(
                    v[h*nfield + i*nlon + j-1],
                    v[h*nfield + i*nlon + j  ],
                    v[h*nfield + i*nlon + j+1],
                    dlambda
                );
            }
            // Eastern boundary
            curl[h*nfield + i*nlon + nlon-1] = curlcoef[i] * fd_kernel_1_eq32(
                v[h*nfield + i*nlon + nlon-2],
                v[h*nfield + i*nlon + nlon-1],
                v[h*nfield + i*nlon         ], // wrapped around periodic boundary
                dlambda
            );
        }
        // Relative vorticity, -du/dphi term
        // Interior domain
        for (size_t i = 1; i < nlat-1; i++) {
            for (size_t j = 0; j < nlon; j++) {
                curl[h*nfield + i*nlon + j] -= curlcoef[i] * fd_kernel_1_eq32(
                    coslat[i-1] * u[h*nfield + (i-1)*nlon + j],
                    coslat[i  ] * u[h*nfield + (i  )*nlon + j],
                    coslat[i+1] * u[h*nfield + (i+1)*nlon + j],
                    dphi
                );
            }
        }
        // Northern boundary if it is the North Pole: replace all values by the
        // mean of the next latitude south. This ensures that there is only one
        // unique value at the pole and avoids div/0 in the curl coefficient.
        if (ispole(lat[0])) {
            const double next_mean = nanmean(curl + (h*nfield + nlon), nlon);
            for (size_t j = 0; j < nlon; j++) {
                ASSIGN_PRESERVE_NAN(curl[h*nfield + j], next_mean);
            }
        // Northern boundary otherwise: use one-sided difference
        } else {
            for (size_t j = 0; j < nlon; j++) {
                curl[h*nfield + j] -= curlcoef[0] * fd_kernel_1_eq32(
                    NAN, // force one-sided difference
                    coslat[0] * u[h*nfield        + j],
                    coslat[1] * u[h*nfield + nlon + j],
                    dphi
                );
            }
        }
        // Southern boundary if it is the South Pole: replace all values by the
        // mean of the next latitude north. This ensures that there is only one
        // unique value at the pole and avoids div/0 in the curl coefficient.
        if (ispole(lat[nlat-1])) {
            const double next_mean = nanmean(curl + (h*nfield + nfield - 2*nlon), nlon);
            for (size_t j = nfield - nlon; j < nfield; j++) {
                ASSIGN_PRESERVE_NAN(curl[h*nfield + j], next_mean);
            }
        // Southern boundary otherwise: use one-sided difference
        } else {
            for (size_t j = 0; j < nlon; j++) {
                curl[h*nfield + (nlat-1)*nlon + j] -= curlcoef[nlat-1] * fd_kernel_1_eq32(
                    coslat[nlat-2] * u[h*nfield + (nlat-2)*nlon + j],
                    coslat[nlat-1] * u[h*nfield + (nlat-1)*nlon + j],
                    NAN, // force one-sided difference
                    dphi
                );
            }
        }
    }
    free(coslat);
    free(curlcoef);
    return 0;
}


void norm_grad_log_abs(
    const double * lat, // [nlat]
    const double * pv, // [nlev, nlat, nlon]
    const double threshold, // TODO
    const size_t nlev,
    const size_t nlat,
    const size_t nlon,
    double * out // [nlev, nlat, nlon]
) {
    const size_t nall = nlev * nlat * nlon;
    for (size_t i = 0; i < nall; i++) {
        const double q = fabs(pv[i]);
        out[i] = (q >= threshold) ? log(q) : NAN;
    }
    double * gx = malloc(nall * sizeof(double));
    double * gy = malloc(nall * sizeof(double));
    horizontal_gradient(lat, out, nlev, nlat, nlon, gx, gy);
    for (size_t i = 0; i < nall; i++) {
        out[i] = sqrt(gx[i]*gx[i] + gy[i]*gy[i]);
    }
    free(gx);
    free(gy);
}


// Potential temperature theta = T * (p0 / p)^kappa
int potential_temperature_isob(
    const double * isob, // [nisob] in: pressure levels
    const double * t, // [nisob, nlat, nlon] in: temperature
    const size_t nisob,
    const size_t nlat,
    const size_t nlon,
    double * th // [nisob, nlat, nlon] out: potential temperature
) {
    const size_t nfield = nlat * nlon;
    for (size_t h = 0; h < nisob; h++) {
        const double factor = pow(P_REF / isob[h], KAPPA);
        for (size_t i = 0; i < nfield; i++) {
            th[h*nfield + i] = t[h*nfield + i] * factor;
        }
    }
    return 0;
}


// Isentropic density sigma = - 1/g dp/dtheta, computed via 1/(dtheta/dp)
int isentropic_density_isob(
    const double * isob, // [nisob] in: pressure levels
    const double * th, // [nisob, nlat, nlon] in: potential temperature
    const size_t nisob,
    const size_t nlat,
    const size_t nlon,
    double * sg // [nisob, nlat, nlon] out: isentropic density
) {
    if (nisob < 3) {
        return 1;
    }
    const size_t nfield = nlat * nlon;
    const double factor = -1. / GRAVITY;
    // Lower boundary: one-sided 2nd order finite difference
    {
        const double dp01 = (isob[1] - isob[0]) * 100.; // hPa -> Pa
        const double dp12 = (isob[2] - isob[1]) * 100.; // hPa -> Pa
        const double coef0 = -(dp01 + dp01 + dp12) / dp01 / (dp01 + dp12);
        const double coef1 = (dp01 + dp12) / (dp01 * dp12);
        const double coef2 = - dp01 / dp12 / (dp01 + dp12);
        for (size_t i = 0; i < nfield; i++) {
            sg[i] = factor / (  coef0 * th[           i]
                              + coef1 * th[  nfield + i]
                              + coef2 * th[2*nfield + i]);
        }
    }
    // Interior domain: centered finite difference
    for (size_t h = 1; h < nisob-1; h++) {
        const double dp01 = (isob[h] - isob[h-1]) * 100.; // hPa -> Pa
        const double dp12 = (isob[h+1] - isob[h]) * 100.; // hPa -> Pa
        const double coef0 = -dp12 / dp01 / (dp01 + dp12);
        const double coef1 = (dp12 - dp01) / (dp01 * dp12);
        const double coef2 = dp01 / dp12 / (dp01 + dp12);
        for (size_t i = 0; i < nfield; i++) {
            sg[h*nfield + i] = factor / (  coef0 * th[(h-1)*nfield + i]
                                         + coef1 * th[ h   *nfield + i]
                                         + coef2 * th[(h+1)*nfield + i]);
        }
    }
    // Upper boundary
    {
        const double dp01 = (isob[nisob-2] - isob[nisob-3]) * 100.; // hPa -> Pa
        const double dp12 = (isob[nisob-1] - isob[nisob-2]) * 100.; // hPa -> Pa
        const double coef0 = dp12 / dp01 / (dp01 + dp12);
        const double coef1 = -(dp01 + dp12) / (dp01 * dp12);
        const double coef2 = (dp12 + dp12 + dp01) / dp12 / (dp01 + dp12);
        for (size_t i = 0; i < nfield; i++) {
            sg[(nisob-1)*nfield + i] = factor / (  coef0 * th[(nisob-3)*nfield + i]
                                                 + coef1 * th[(nisob-2)*nfield + i]
                                                 + coef2 * th[(nisob-1)*nfield + i]);
        }
    }
    // Set density for gridpoints with unstable stratification to zero
    for (size_t i = 0; i < nisob*nlat*nlon; i++) {
        if (isunstable(sg[i])) {
            sg[i] = 0.;
        }
    }
    return 0;
}


// Mask regions where data has been extrapolated into the ground
int mask_underground_isob(
    const double * isob, // [nisob] in: pressure levels
    const double * psfc, // [nlat, nlon] in: surface pressure
    double * field, // [nisob, nlat, nlon] in/out
    const double fill,
    const size_t nisob,
    const size_t nlat,
    const size_t nlon
) {
    const size_t nfield = nlat * nlon;
    for (size_t h = 0; h < nisob; h++) {
        for (size_t i = 0; i < nfield; i++) {
            if (isob[h] > psfc[i]) {
                field[h*nfield+i] = fill;
            }
        }
    }
    return 0;
}


// Interpolation from pressure to isentropic levels (isentropic density and
// horizontal wind components)
int interpolate_isob_to_isen(
    const double * isob, // [nisob] in: pressure levels (ascending, top-down)
    const double * isen, // [nisen] in: isentropic levels (descending, top-down)
    const double * th, // [nisob, nlat, nlon] in: potential temperature
    const size_t nisob,
    const size_t nisen,
    const size_t nlat,
    const size_t nlon,
    const size_t nvar,
    double * p, // [nisen, nlat, nlon] out: interpolated pressure
    ... // [nisob, nlat, nlon], [nisen, nlat, nlon], ... alternating in/out
) {
    if (nisob < 2) {
        return 1;
    }
    const size_t nfield = nlat * nlon;
    // Take the alternating input/output fields from the variadic arguments and
    // initialize the output fields with NaN values
    double * var[nvar];
    double * var_out[nvar];
    va_list args;
    va_start(args, p);
    for (size_t a = 0; a < nvar; a++) {
        var[a]  = va_arg(args, double *);
        var_out[a] = va_arg(args, double *);
        for (size_t i = 0; i < nisen*nfield; i++) {
            var_out[a][i] = NAN;
        }
    }
    va_end(args);
    // Also initialize the output pressure field with NaN
    for (size_t i = 0; i < nisen*nfield; i++) {
        p[i] = NAN;
    }
    // Interpolate column-wise at every gridpoint
    for (size_t i = 0; i < nfield; i++) {
        size_t h = 0; // index of current pressure level
        size_t m = 0; // index of current isentropic level
        while (h < nisob-1 && m < nisen) {
            const double th0 = th[ h   *nfield + i];
            const double th1 = th[(h+1)*nfield + i];
            // Isentrope in the current pressure interval: interpolate, go to
            // next isentrope, stay in current pressure interval
            if (isin(isen[m], th0, th1)) {
                // Interpolation coefficients (linear)
                const double coeff0 = (th1 - isen[m]) / (th1 - th0);
                const double coeff1 = (isen[m] - th0) / (th1 - th0);
                // Interpolate pressure
                p[m*nfield + i] = coeff0 * isob[h] + coeff1 * isob[h+1];
                // Interpolate other fields
                for (size_t a = 0; a < nvar; a++) {
                    var_out[a][m*nfield + i] =   coeff0 * var[a][ h   *nfield + i]
                                               + coeff1 * var[a][(h+1)*nfield + i];
                }
                m++;
            // Isentrope above the current pressure interval: out-of-range
            // isentrope (top-down scan), keep NaN values, go to next isentrope
            } else if (th0 < isen[m] && th1 < isen[m]) {
                m++;
            // Isentrope below the current pressure interval: go to next
            // pressure interval, check isentrope again
            } else {
                h++;
            }
        }
        // Remaining isentropes are out-of-range: keep NaN values
    }
    return 0;
}


// Mask regions where data has been extrapolated into the ground (in-place)
int mask_underground_isen(
    const double * p, // [nisen, nlat, nlon] in: pressure on isentropes
    const double * psfc, // [nlat, nlon] in: surface pressure
    double * field, // [nisen, nlat, nlon] in/out
    const double fill,
    const size_t nisen,
    const size_t nlat,
    const size_t nlon
) {
    size_t nfield = nlat * nlon;
    for (size_t m = 0; m < nisen; m++) {
        for (size_t i = 0; i < nfield; i++) {
            if (p[m*nfield + i] > psfc[i]) {
                field[m*nfield + i] = fill;
            }
        }
    }
    return 0;
}


// Absolute vorticity calculation is vertical coordinate-agnostic
int absolute_vorticity(
    const double * lat, // [nlat] in: latitude grid
    const double * u, // [nlev, nlat, nlon] in: zonal wind
    const double * v, // [nlev, nlat, nlon] in: meridional wind
    const size_t nlev,
    const size_t nlat,
    const size_t nlon,
    double * av // [nlev, nlat, nlon] out: absolute vorticity
) {
    size_t nfield = nlat * nlon;
    // Compute relative vorticity
    int status = curl(lat, u, v, nlev, nlat, nlon, av);
    if (status != 0) {
        return 10 + status;
    }
    // Planetary vorticity
    double * coriolis = malloc(nlat * sizeof(double)); // planetary vorticity
    for (size_t i = 0; i < nlat; ++i) {
        coriolis[i] = 2. * OMEGA * sin(deg2rad(lat[i]));
    }
    // Add the planetary vorticity term. If wind or relative vorticity is
    // undefined (NaN) in a cell, put the planetary vorticity anyway.
    for (size_t h = 0; h < nlev; h++) {
        for (size_t i = 0; i < nlat; i++) {
            for (size_t j = 0; j < nlon; j++) {
                const size_t idx = h*nfield + i*nlon + j;
                if (isnan(av[idx])) {
                    av[idx] = coriolis[i];
                } else {
                    av[idx] += coriolis[i];
                }
            }
        }
    }
    free(coriolis);
    return 0;
}


int potential_vorticity_isen(
    const double * av, // [nisen, nlat, nlon] in: absolute vorticity
    const double * sg, // [nisen, nlat, nlon] in: isentropic density
    const size_t nisen,
    const size_t nlat,
    const size_t nlon,
    double * pv // [nisen, nlat, nlon] out: potential vorticity
) {
    const size_t ndomain = nisen * nlat * nlon;
    for (size_t i = 0; i < ndomain; i++) {
        // PV is set to NaN for underground and unstable gridpoints
        pv[i] = isstable(sg[i]) ? (av[i] / sg[i] * PVU) : NAN;
    }
    return 0;
}


void isob_to_isen_all(
    const double * isob, // [nisob] in: pressure levels
    const double * isen, // [nisen] in: isentropes
    const double * lat, // [nlat] in: latitude grid
    const double * t_isob, // [nvec, nisob, nlat, nlon] in: temperature
    const double * u_isob, // [nvec, nisob, nlat, nlon] in: zonal wind
    const double * v_isob, // [nvec, nisob, nlat, nlon] in: meridional wind
    const double * psfc, // NULL or [nvec, nlat, nlon] in: surface pressure
    const size_t nvec,
    const size_t nisob,
    const size_t nisen,
    const size_t nlat,
    const size_t nlon,
    double * p_isen, // [nvec, nisen, nlat, nlon] out: pressure
    double * u_isen, // [nvec, nisen, nlat, nlon] out: zonal wind
    double * v_isen, // [nvec, nisen, nlat, nlon] out: meridional wind
    double * sg_isen, // [nvec, nisen, nlat, nlon] out: isentropic density
    double * av_isen, // [nvec, nisen, nlat, nlon] out: absolute vorticity
    double * pv_isen // [nvec, nisen, nlat, nlon] out: potential vorticity
) {
    const size_t nfld = nlat * nlon;
    const size_t nchunk_isob = nisob * nlat * nlon;
    const size_t nchunk_isen = nisen * nlat * nlon;
    for (size_t k = 0; k < nvec; k++) {
        // Shorthands for current input fields
        const double * tt_isob = t_isob + (k * nchunk_isob);
        const double * uu_isob = u_isob + (k * nchunk_isob);
        const double * vv_isob = v_isob + (k * nchunk_isob);
        //const double * ps_isob = psfc + (k * nfld);
        // Shorthands for current output fields
        double * p  = p_isen  + (k * nchunk_isen);
        double * u  = u_isen  + (k * nchunk_isen);
        double * v  = v_isen  + (k * nchunk_isen);
        double * sg = sg_isen + (k * nchunk_isen);
        double * av = av_isen + (k * nchunk_isen);
        double * pv = pv_isen + (k * nchunk_isen);
        // Temporary variables for pressure-level computations
        double * th_isob = malloc(nchunk_isob * sizeof(double));
        double * sg_isob = malloc(nchunk_isob * sizeof(double));
        // Potential temperature θ = θ(p, T)
        potential_temperature_isob(isob, tt_isob, nisob, nlat, nlon, th_isob);
        // Isentropic density σ = σ(p, θ)
        isentropic_density_isob(isob, th_isob, nisob, nlat, nlon, sg_isob);
        // Interpolate p, u, v and σ to isentropes
        interpolate_isob_to_isen(
            isob, isen, th_isob,
            nisob, nisen, nlat, nlon, 3, // three fields other then pressure: u, v, sg
            p, uu_isob, u, vv_isob, v, sg_isob, sg
        );
        // Everything is now in isentropic coordinates, free isobaric variables
        free(th_isob);
        free(sg_isob);
        // If a surface pressure field is available, fill u and v with NaN and
        // σ with 0 underground
        if (psfc) {
            const double * ps = psfc + k * nfld;
            mask_underground_isen(p, ps, u, NAN, nisen, nlat, nlon);
            mask_underground_isen(p, ps, v, NAN, nisen, nlat, nlon);
            mask_underground_isen(p, ps, sg, 0., nisen, nlat, nlon);
        }
        // Absolute vorticity ω = ω(u, v)
        absolute_vorticity(lat, u, v, nisen, nlat, nlon, av);
        // Potential vorticity q = q(ω, σ)
        potential_vorticity_isen(av, sg, nisen, nlat, nlon, pv);
    }
}

