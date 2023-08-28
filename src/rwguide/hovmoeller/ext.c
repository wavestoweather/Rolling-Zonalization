#include <math.h>
#include <stdlib.h>
#include <stdbool.h>


// Is x contained in the interval [x0, x1) if x0 < x1 or [x1, x0) if x1 < x0?
inline bool isin(const double x, const double x0, const double x1) {
    return (x0 <= x) ? (x < x1) : (x1 <= x);
}

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
) {
    const size_t nfld = nlat * nlon;
    // Split latitude kernel in two for centered application
    const size_t nlkn = nlat_kernel / 2;
    const size_t nlks = nlat_kernel - nlkn;
    // Levels are independent
    for (size_t h = 0; h < nlev; h++) {
        const double * x = field_x + h*nfld;
        const double * y = field_y + h*nfld;
        // Weighted mean aggregators
        double * ys = malloc((nlon + nlon_kernel) * sizeof(double)); // extend for periodicity
        double * ws = malloc((nlon + nlon_kernel) * sizeof(double)); // extend for periodicity
        // Aggregation along meridians
        for (size_t j = 0; j < nlon; j++) {
            size_t idx = nlat; // fallback if contour is not found
            // Find index of contour
            for (size_t i = 0; i < nlat-1; i++) {
                if (isin(xvalue, x[i*nlon + j], x[(i+1)*nlon + j])) {
                    idx = i;
                    break;
                }
            }
            // Contour not found at this longitude, put NaN and move on
            if (idx == nlat) {
                ys[j] = NAN;
                ws[j] = NAN;
                continue;
            }
            // Initialize weighted mean aggregators
            ys[j] = 0.; // sum of weighted values
            ws[j] = 0.; // sum of weights
            // Meridional extent of averaging
            size_t idxn = (idx < nlkn) ? 0 : idx - nlkn; // clip at 0
            size_t idxs = ((idx + nlks) > (nlat - 1)) ? nlat-1 : idx + nlks; // clip at nlat-1
            for (size_t i = idxn; i < idxs; i++) {
                const double weight = lat_kernel[i-idxn] * area_weights[i];
                ys[j] += y[i*nlon + j] * weight;
                ws[j] += weight;
            }
        }
        // Repeat values to simplify handling of periodicity
        for (size_t j = 0; j < nlon_kernel; j++) {
            ys[nlon + j] = ys[j];
            ws[nlon + j] = ws[j];
        }
        // Compute weighted sum in lon (with periodicity)
        for (size_t j = 0; j < nlon; j++) {
            const size_t center = (j + nlon_kernel / 2) % nlon;
            double sum = 0.;
            double wsum = 0.;
            for (size_t k = 0; k < nlon_kernel; k++) {
                if (!isnan(ys[j+k])) {
                    sum += ys[j+k] * lon_kernel[k];
                    wsum += ws[j+k] * lon_kernel[k];
                }
            }
            mean_y[h*nlon + center] = (wsum > 0) ? sum / wsum : NAN;
        }
    }
}

