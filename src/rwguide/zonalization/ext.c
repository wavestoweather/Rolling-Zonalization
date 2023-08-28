#include <math.h>
#include <stdlib.h>
#include <stddef.h> 
#include <stdbool.h>

// Smallest isentropic density value still considered stable and above ground
#define PVU 1.0e6 // multiply to obtain PV in PVU
#define TOL_SG 0.001 // kg / K / m²

// Is the gridbox stable based on isentropic density?
inline bool isstable(const double sg) {
    return TOL_SG <= sg;
}

inline void fill_nan(double * arr, const size_t n) {
    for (size_t i = 1; i < n; i++) {
        arr[i] = NAN;
    }
}

inline void cumulative_sum(double * arr, const size_t n) {
    for (size_t i = 1; i < n; i++) {
        arr[i] += arr[i-1];
    }
}

inline void linspace(const double start, const double end, double * arr, const size_t n) {
    for (size_t k = 0; k < n; ++k) {
        arr[k] = start + (end - start) * ((double) k) / ((double) (n - 1));
    }
}

inline void interp_desc(
    const double * xs, // where to interpolate to
    const double * xp, // coordinates...
    const double * fp, // ...and corresponding values
    const double fp_max, // value to use for max overflow
    const double fp_min, // value to use for min overflow
    const size_t ns,
    const size_t np,
    double * fs // output: interpolated function
) {
    size_t s = 0;
    size_t p = 0;
    // Fill until first value in interval
    while (xs[s] < xp[p]) {
        fs[s] = fp_max;
        s++;
    }
    while (p < (np+1) && s < ns) {
        // Go to next interval if value out of range or interval is empty
        if (xs[s] > xp[p+1] || xp[p] == xp[p+1]) {
            p++;
        // Linear interpolation in interval
        } else {
            const double dx = xp[p+1] - xp[p];
            const double dy = fp[p+1] - fp[p];
            fs[s] = fp[p] + (xs[s] - xp[p]) * dy / dx;
            s++;
        }
    }
    // Fill until end
    while (s < ns) {
        fs[s] = fp_min;
        s++;
    }
}


/* Background density field computation */

inline double weighted_nanmean_1d(
    const double * weights,
    const double * signal,
    const double threshold,
    const size_t n
) {
    double signal_agg = 0.;
    double weight_agg = 0.;
    for (size_t i = 0; i < n; i++) {
        if (!isnan(signal[i])) {
            signal_agg += signal[i] * weights[i];
            weight_agg += weights[i]; // sum weights for normalization
        }
    }
    // If the weights aggregated for all valid gridpoints do not sum to at
    // least the given threshold value, the result is not considered to be
    // representative.
    if (weight_agg <= threshold) {
        return NAN;
    }
    // Divide by the sum of all weights applied to valid values. This probably
    // does not correspond to an ideal redistribution of weights, but it's
    // easy. It would be better to redistribute the unused weights
    // preferentially to nearby gridpoints, instead of uniformly.
    return signal_agg / weight_agg;
}

int weighted_nanmean_rolling_1d(
    const double * weights,
    const double * signal,
    const double threshold,
    const size_t nsignal,
    const size_t nweights,
    double * out // same size as signal
) {
    // Windows should be folded first
    if (nweights > nsignal) {
        return 1;
    }
    // Create a wrapped copy of the signal
    double * wrapped_signal = malloc((nsignal + nweights) * sizeof(double));
    for (size_t i = 0; i < nsignal; i++) {
        wrapped_signal[i] = signal[i];
    }
    for (size_t i = 0; i < nweights; i++) {
        wrapped_signal[nsignal + i] = signal[i];
    }
    // Compute the weighted rolling mean
    for (size_t i = 0; i < nsignal; i++) {
        const size_t center = (i + nweights / 2) % nsignal;
        out[center] = weighted_nanmean_1d(
            weights,
            wrapped_signal + i, // shift to starting position
            threshold,
            nweights
        );
    }
    free(wrapped_signal);
    return 0;
}


int rolling_mean_background(
    const double * window,
    const double * field,
    const size_t nlev,
    const size_t nlat,
    const size_t nlon,
    const size_t nwin,
    double * out
) {
    const size_t nfield = nlat * nlon;
    for (size_t h = 0; h < nlev; h++) {
        for (size_t i = 0; i < nlat; i++) {
            int status = weighted_nanmean_rolling_1d(
                window + (i * nwin),
                field + (h * nfield) + (i * nlon),
                0., // threshold
                nlon,
                nwin,
                out + (h * nfield) + (i * nlon)
            );
            if (status != 0) {
                return status;
            }
        }
    }
    return 0;
}


/* Quadrature */

void latitude_mass_north_2d(
    const double * area, // [nlat, nlon]
    const double * sg, // [nlat, nadv]
    const size_t nlat,
    const size_t nlon,
    const size_t nadv,
    double * mass // [nlat]
) {
    // Initialize all mass accumulators with zero
    for (size_t i = 0; i < nlat; i++) {
        mass[i] = 0.;
    }
    // No mass north of northernmost latitude, keep zero value. The area
    // weights for the northernmost latitude should already account for the
    // fact that the gridbox only exists south of the associated latitude. No
    // need to split the gridbox further here.
    for (size_t j = 0; j < nlon; j++) {
        mass[1] += sg[j] * area[j];
    }
    // The area weights for the inner domain should represent the entire
    // associated gridbox. Split them in half here and attribute the northern
    // half to the mass associated with the current latitude and the southern
    // half to the mass associated with the next latitude.
    for (size_t i = 1; i < nlat-1; i++) {
        for (size_t j = 0; j < nlon; j++) {
            // Account for periodicity-duplicated values with nadv and split
            // the gridbox
            mass[i  ] += 0.5 * sg[i*nadv + j] * area[i*nlon + j];
            mass[i+1] += 0.5 * sg[i*nadv + j] * area[i*nlon + j];
        }
    }
    // The area weights for the southernmost latitude should already account
    // for the fact that the gridbox only exists north of the associated
    // latitude. No need to split the gridbox further here.
    for (size_t j = 0; j < nlon; j++) {
        mass[nlat-1] += sg[(nlat-1)*nadv + j] * area[(nlat-1)*nlon + j];
    }
    // So far only mass between latitudes computed, now aggregate
    cumulative_sum(mass, nlat);
}


void contour_mass_north_2d(
    const double * area, // [nlat, nlon]
    const double * av, // [nlat, nadv]
    const double * sg, // [nlat, nadv]
    const double * pvc, // [npvc] descending order
    const size_t nlat,
    const size_t nlon,
    const size_t nadv,
    const size_t npvc,
    double * mass // [npvc]
) {
    // Initialize all mass accumulators with zero
    for (size_t k = 0; k < npvc; k++) {
        mass[k] = 0.;
    }
    for (size_t i = 0; i < nlat; i++) {
        for (size_t j = 0; j < nlon; j++) {
            const double weight = area[i*nlon + j];
            if (weight == 0.) {
                continue;
            }
            // Only consider gridboxes with stable stratification
            const double density = sg[i*nadv + j];
            if (isstable(density)) {
                // Ertel potential vorticity
                const double pv = av[i*nadv + j] / density * PVU;
                // Compute mass between contours, add mass of gridbox to
                // accumulator for closest contour
                for (size_t k = 0; k < npvc; k++) {
                    if (pvc[k] <= pv) {
                        mass[k] += density * weight;
                        break;
                    }
                }
                // If no match found above, skip the gridbox
            }
        }
    }
    // So far only mass between contours computed, now aggregate
    cumulative_sum(mass, npvc);
}


/* Zonalization */

void zonalize_north_2d(
    const double * area, // [nlat, nlon] area weights for quadrature
    const double * av, // [nlat, nadv] absolute vorticity
    const double * sg, // [nlat, nadv] isentropic density
    const double * bg, // [nlat, nadv] isentropic density (background)
    const size_t nlat,
    const size_t nlon,
    const size_t nadv,
    double * pv_zonalized // [nlat]
) {
    // Determine area north of latitudes of background: use background
    // isentropic density field
    double * latitude_mass = malloc(nlat * sizeof(double));
    latitude_mass_north_2d(area, bg, nlat, nlon, nadv, latitude_mass);
    // Find maximum of PV for contour generation, use 0 for min.
    // Assumption: northern hemisphere.
    double max_pv = 0.;
    for (size_t i = 0; i < nlat; i++) {
        for (size_t j = 0; j < nlon; j++) {
            const double density = sg[i*nadv + j];
            if (isstable(density)) {
                const double pv = av[i*nadv + j] / density * PVU;
                if (pv > max_pv) {
                    max_pv = pv;
                }
            }
         }
    }
    // If no valid PV values were found, fill output with NaN and return
    if (max_pv == 0.) {
        fill_nan(pv_zonalized, nlat);
        return;
    }
    // Generate 10 contours (desc) to estimate the zonalized PV profile:
    // determine area north of contours for equivalent latitude
    const size_t npvc_approx = 10;
    double * pvc_approx          = malloc(npvc_approx * sizeof(double));
    double * contour_mass_approx = malloc(npvc_approx * sizeof(double));
    linspace(max_pv, 0., pvc_approx, npvc_approx);
    contour_mass_north_2d(area, av, sg, pvc_approx, nlat, nlon, nadv, npvc_approx, contour_mass_approx);
    // Estimate zonalized PV on the input latitude grid. These contour values
    // are used as a refined contour selection in a second pass.
    double * pvc_refine = malloc(nlat * sizeof(double));
    interp_desc(latitude_mass, contour_mass_approx, pvc_approx, max_pv, 0., nlat, npvc_approx, pvc_refine);
    free(pvc_approx);
    free(contour_mass_approx);
    // The refined contour set is filled up with 0 after the interpolation
    // range, figure out where the valid contours end for the second pass
    // (avoid multiple 0-contours)
    size_t npvc_refine = 0;
    while (npvc_refine < nlat && pvc_refine[npvc_refine] > 0.) {
        npvc_refine++;
    }
    // Determine area north of refined contours for equivalent latitude
    double * contour_mass_refine = malloc(npvc_refine * sizeof(double));
    contour_mass_north_2d(area, av, sg, pvc_refine, nlat, nlon, nadv, npvc_refine, contour_mass_refine);
    // Interpolate PV onto latitudes via area, fill with NaN if out of range
    interp_desc(latitude_mass, contour_mass_refine, pvc_refine, max_pv, NAN, nlat, npvc_refine, pv_zonalized);
    free(pvc_refine);
    free(contour_mass_refine);
    free(latitude_mass);
}


// Zonalization for data with descending latitude
void zonalize(
    const double * area, // [nlat, nlon] area weights for quadrature
    const double * av, // [nlev, nlat, nlon] absolute vorticity
    const double * sg, // [nlev, nlat, nlon] isentropic density
    const double * bg, // [nlev, nlat, nlon] isentropic density (background)
    const size_t nh_last, // index of northern hemisphere end (NH in [0, nh_last])
    const size_t sh_first, // index of southern hemisphere start (SH in [sh_first, nlat-1])
    const size_t nlev,
    const size_t nlat,
    const size_t nlon,
    double * pv_zonalized // [nlev, nlat]
) {
    // Allow specification of NH only data with inh >= nlat (crop at nlat-1)
    const size_t inh = (nh_last > nlat-1) ? nlat-1 : nh_last;
    // Allow specification of SH only data with ish <= 0 (crop at 0)
    const size_t ish = (sh_first < 0) ? 0 : sh_first;
    // Determine sizes of northern and southern hemisphere
    const size_t nnh = inh + 1;
    const size_t nsh = nlat - ish;
    // Create a latitude-flipped copy of area for the southern hemisphere
    double * area_sh = malloc((nsh * nlon) * sizeof(double));
    for (size_t i = 0; i < nsh; i++) {
        for (size_t j = 0; j < nlon; j++) {
            const size_t iflip = nsh - i - 1;
            area_sh[iflip*nlon + j] = area[(ish+i)*nlon + j];
        }
    }
    // Zonalize each level independently
    const size_t nfld = nlat * nlon;
    for (size_t h = 0; h < nlev; h++) {
        // Northern Hemisphere
        if (nnh > 1) {
            const double * av_nh = av + h*nfld;
            const double * sg_nh = sg + h*nfld;
            const double * bg_nh = bg + h*nfld;
            double * pvz_nh = pv_zonalized + h*nlat;
            zonalize_north_2d(area, av_nh, sg_nh, bg_nh, nnh, nlon, nlon, pvz_nh);
        }
        // Fill values between hemispheres with NaN
        for (size_t i = inh+1; i < ish; i++) {
            pv_zonalized[h*nlat + i] = NAN;
        }
        // Southern Hemisphere: flip fields on the equator and multiply
        // absolute vorticity by -1 to obtain a NH field, then use the NH
        // zonalization and return the results back to SH values and ordering
        if (nsh > 1) {
            // Create a copy with flipped latitude axis so data starts with
            // south pole
            double * av_sh = malloc((nsh * nlon) * sizeof(double));
            double * sg_sh = malloc((nsh * nlon) * sizeof(double));
            double * bg_sh = malloc((nsh * nlon) * sizeof(double));
            for (size_t i = 0; i < nsh; i++) {
                const size_t iflip = nsh - i - 1;
                for (size_t j = 0; j < nlon; j++) {
                    // Flip sign of absolute vorticity to obtain NH-like values
                    av_sh[iflip*nlon + j] = -av[h*nfld + (ish+i)*nlon + j];
                    sg_sh[iflip*nlon + j] =  sg[h*nfld + (ish+i)*nlon + j];
                    bg_sh[iflip*nlon + j] =  bg[h*nfld + (ish+i)*nlon + j];
                }
            }
            double * temp = malloc(nsh * sizeof(double));
            // Periodicity-duplicated data is set up to resemble a northern
            // hemisphere state, so NH zonalization can be used
            zonalize_north_2d(area_sh, av_sh, sg_sh, bg_sh, nsh, nlon, nlon, temp);
            for (size_t i = 0; i < nsh; i++) {
                const size_t iflip = nsh - i - 1;
                // Flip sign and latitude back when assigning zonalized PV
                pv_zonalized[h*nlat + ish + i] = -temp[iflip];
            }
            free(temp);
            free(av_sh);
            free(sg_sh);
            free(bg_sh);
        }
    }
    free(area_sh);
}


// Rolling-window zonalization for data with descending latitude
void zonalize_rolling(
    const double * area, // [nlat, nwin] area weights for quadrature
    const double * av, // [nlev, nlat, nlon] absolute vorticity
    const double * sg, // [nlev, nlat, nlon] isentropic density
    const double * bg, // [nlev, nlat, nlon] isentropic density (background)
    const size_t nh_last, // index of northern hemisphere end (NH in [0, nh_last])
    const size_t sh_first, // index of southern hemisphere start (SH in [sh_first, nlat-1])
    const size_t nlev,
    const size_t nlat,
    const size_t nlon,
    const size_t nwin,
    double * pv_zonalized // [nlev, nlat, nlon]
) {
    // Allow specification of NH only data with inh >= nlat (crop at nlat-1)
    const size_t inh = (nh_last > nlat-1) ? nlat-1 : nh_last;
    // Allow specification of SH only data with ish <= 0 (crop at 0)
    const size_t ish = (sh_first < 0) ? 0 : sh_first;
    // Determine sizes of northern and southern hemisphere
    const size_t nnh = inh + 1;
    const size_t nsh = nlat - ish;
    // Create a latitude-flipped copy of area for the southern hemisphere
    double * area_sh = malloc((nsh * nwin) * sizeof(double));
    for (size_t i = 0; i < nsh; i++) {
        for (size_t j = 0; j < nwin; j++) {
            const size_t iflip = nsh - i - 1;
            area_sh[iflip*nwin + j] = area[(ish+i)*nwin + j];
        }
    }
    // Zonalize each level independently
    const size_t nfld = nlat * nlon;
    const size_t nadv = nlon + nwin; // for periodicity-duplicated domain
    for (size_t h = 0; h < nlev; h++) {
        // Northern Hemisphere
        if (nnh > 1) {
            // Create a copy with duplicated values for periodicity
            double * av_nh = malloc((nnh * nadv) * sizeof(double));
            double * sg_nh = malloc((nnh * nadv) * sizeof(double));
            double * bg_nh = malloc((nnh * nadv) * sizeof(double));
            for (size_t i = 0; i < nnh; i++) {
                for (size_t j = 0; j < nlon; j++) {
                    av_nh[i*nadv + j] = av[h*nfld + i*nlon + j];
                    sg_nh[i*nadv + j] = sg[h*nfld + i*nlon + j];
                    bg_nh[i*nadv + j] = bg[h*nfld + i*nlon + j];
                }
                for (size_t j = 0; j < nwin; j++) {
                    av_nh[i*nadv + nlon + j] = av[h*nfld + i*nlon + j];
                    sg_nh[i*nadv + nlon + j] = sg[h*nfld + i*nlon + j];
                    bg_nh[i*nadv + nlon + j] = bg[h*nfld + i*nlon + j];
                }
            }
            for (size_t j = 0; j < nlon; j++) {
                double * temp = malloc(nnh * sizeof(double));
                zonalize_north_2d(area, av_nh+j, sg_nh+j, bg_nh+j, nnh, nwin, nadv, temp);
                const size_t center = (j + nwin / 2) % nlon;
                for (size_t i = 0; i < nnh; i++) {
                    pv_zonalized[h*nfld + i*nlon + center] = temp[i];
                }
                free(temp);
            }
            free(av_nh);
            free(sg_nh);
            free(bg_nh);
        }
        // Fill values between hemispheres with NaN
        for (size_t i = inh+1; i < ish; i++) {
            for (size_t j = 0; j < nlon; j++) {
                pv_zonalized[h*nfld + i*nlon + j] = NAN;
            }
        }
        // Southern Hemisphere: flip fields on the equator and multiply
        // absolute vorticity by -1 to obtain a NH field, then use the NH
        // zonalization and return the results back to SH values and ordering
        if (nsh > 1) {
            // Create a copy with duplicated values for periodicity
            double * av_sh = malloc((nsh * nadv) * sizeof(double));
            double * sg_sh = malloc((nsh * nadv) * sizeof(double));
            double * bg_sh = malloc((nsh * nadv) * sizeof(double));
            for (size_t i = 0; i < nsh; i++) {
                // Flip latitude axis so data starts with pole
                const size_t iflip = nsh - i - 1;
                for (size_t j = 0; j < nlon; j++) {
                    // Flip sign of absolute vorticity to obtain NH-like values
                    av_sh[iflip*nadv + j] = -av[h*nfld + (ish+i)*nlon + j];
                    sg_sh[iflip*nadv + j] =  sg[h*nfld + (ish+i)*nlon + j];
                    bg_sh[iflip*nadv + j] =  bg[h*nfld + (ish+i)*nlon + j];
                }
                for (size_t j = 0; j < nwin; j++) {
                    av_sh[iflip*nadv + nlon + j] = -av[h*nfld + (ish+i)*nlon + j];
                    sg_sh[iflip*nadv + nlon + j] =  sg[h*nfld + (ish+i)*nlon + j];
                    bg_sh[iflip*nadv + nlon + j] =  bg[h*nfld + (ish+i)*nlon + j];
                }
            }
            for (size_t j = 0; j < nlon; j++) {
                double * temp = malloc(nsh * sizeof(double));
                // Periodicity-duplicated data is set up to resemble a northern
                // hemisphere state, so NH zonalization can be used
                zonalize_north_2d(area_sh, av_sh+j, sg_sh+j, bg_sh+j, nsh, nwin, nadv, temp);
                const size_t center = (j + nwin / 2) % nlon;
                for (size_t i = 0; i < nsh; i++) {
                    const size_t iflip = nsh - i - 1;
                    // Flip sign and latitude back when assigning zonalized PV
                    pv_zonalized[h*nfld + (ish+i)*nlon + center] = -temp[iflip];
                }
                free(temp);
            }
            free(av_sh);
            free(sg_sh);
            free(bg_sh);
        }
    }
    free(area_sh);
}

