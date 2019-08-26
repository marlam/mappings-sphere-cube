#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>

#include "spherecube.hpp"

#include "common.hpp"
#include "pngvis.hpp"


template <typename FLOAT>
void analyze(const Context<FLOAT>& ctx, int w, int h, FLOAT* dist_area, FLOAT* dist_isot, FLOAT* snyder_s, FLOAT* snyder_omega)
{
    const FLOAT area_unit_sphere = 4 * m_pi<FLOAT>();
    const FLOAT area_unit_square = 4;
    const FLOAT R = area_unit_square / (area_unit_sphere / 6);
    for (int iy = 0; iy < h; iy++) {
        FLOAT y = ((iy + static_cast<FLOAT>(0.5)) / h) * 2 - 1;
        for (int ix = 0; ix < w; ix++) {
            FLOAT x = ((ix + static_cast<FLOAT>(0.5)) / w) * 2 - 1;
            FLOAT lat, lon;
            ctx.inverse(x, y, lat, lon);
            if (std::isfinite(lat) && std::isfinite(lon)) {
                FLOAT h, k, theta_prime, a, b, omega, s;
                ctx.analysis(lat, lon, h, k, theta_prime, a, b, omega, s);
                dist_area[iy * w + ix] = s / R;
                dist_isot [iy * w + ix] = a / b;
                snyder_s[iy * w + ix] = s;
                snyder_omega[iy * w + ix] = degrees(omega);
            } else {
                dist_area[iy * w + ix] = NAN;
                dist_isot [iy * w + ix] = NAN;
                snyder_s[iy * w + ix] = NAN;
                snyder_omega[iy * w + ix] = NAN;
            }
        }
    }
}

template <typename FLOAT>
void scan_array(const FLOAT* data, int size,
        FLOAT& min, FLOAT& min_1percent,
        FLOAT& max, FLOAT& max_1percent,
        FLOAT& mean, FLOAT& median,
        FLOAT& max_abs_error, FLOAT& rmse)
{
    std::vector<FLOAT> valid_data;
    valid_data.reserve(size);
    min = max = NAN;
    mean = 0;
    rmse = 0;
    for (int i = 0; i < size; i++) {
        if (std::isfinite(data[i])) {
            if (!std::isfinite(min) || data[i] < min)
                min = data[i];
            if (!std::isfinite(max) || data[i] > max)
                max = data[i];
            mean += data[i];
            FLOAT abs_err = std::fabs(1 - data[i]);
            rmse += abs_err * abs_err;
            valid_data.push_back(data[i]);
        }
    }
    if (valid_data.size() == 0) {
        mean = NAN;
        min_1percent = NAN;
        max_1percent = NAN;
        median = NAN;
        max_abs_error = NAN;
        rmse = NAN;
    } else {
        mean /= valid_data.size();
        max_abs_error = std::max(1 - min, max - 1);
        rmse = std::sqrt(rmse / valid_data.size());
        std::sort(valid_data.begin(), valid_data.end());
        if (valid_data.size() % 2 == 0) {
            median = (valid_data[valid_data.size() / 2 - 1]
                    + valid_data[valid_data.size() / 2]) / 2;
        } else {
            median = valid_data[valid_data.size() / 2];
        }
        min_1percent = valid_data[valid_data.size() / 100];
        max_1percent = valid_data[valid_data.size() - 1 - valid_data.size() / 100];
    }
}

/*
 * Main function.
 */

int main(void)
{
    const int w = 512;
    const int h = 512;

    double *dist_area_d = new double[w * h];
    double *dist_isot_d = new double[w * h];
    double *s_d = new double[w * h];
    double *omega_d = new double[w * h];
    double min_d, max_d; //, mean_d;
    double min_1percent_d, max_1percent_d, median_d;

    double da_mean_d, da_max_abs_error_d, da_rmse_d;
    double di_mean_d, di_max_abs_error_d, di_rmse_d;

    // Loop over the projection methods
    for (int method = 0; method < n_projections; method++) {
        Context<double> ctx_d(projections[method].proj,
                projections[method].params,
                projections[method].center);

        fprintf(stderr, "%s double\n", projections[method].name);
        analyze<double>(ctx_d, w, h, dist_area_d, dist_isot_d, s_d, omega_d);

#if 0
        scan_array(s_d, w * h, min_d, min_1percent_d, max_d, max_1percent_d, mean_d, median_d, max_abs_error_d, rmse_d);
        fprintf(stderr, "  Snyder s        min=%g min_1percent=%g max=%g max_1percent=%g mean=%g median=%g\n", min_d, min_1percent_d, max_d, max_1percent_d, mean_d, median_d);
        fprintf(stderr, "  Snyder s        min=%g min_1percent=%g max=%g max_1percent=%g mean=%g median=%g\n", min_d, min_1percent_d, max_d, max_1percent_d, mean_d, median_d);
        scan_array(omega_d, w * h, min_d, min_1percent_d, max_d, max_1percent_d, mean_d, median_d, max_abs_error_d, rmse_d);
        fprintf(stderr, "  Snyder omega    min=%g min_1percent=%g max=%g max_1percent=%g mean=%g median=%g\n", min_d, min_1percent_d, max_d, max_1percent_d, mean_d, median_d);
        fprintf(stderr, "  Snyder omega    min=%g min_1percent=%g max=%g max_1percent=%g mean=%g median=%g\n", min_d, min_1percent_d, max_d, max_1percent_d, mean_d, median_d);
#endif

        scan_array(dist_area_d, w * h, min_d, min_1percent_d, max_d, max_1percent_d, da_mean_d, median_d, da_max_abs_error_d, da_rmse_d);
        scan_array(dist_isot_d, w * h, min_d, min_1percent_d, max_d, max_1percent_d, di_mean_d, median_d, di_max_abs_error_d, di_rmse_d);

        fprintf(stderr, "  area dist (D_A) avg | max abs | rmse : %.2f | %.2f | %.2f\n", da_mean_d, da_max_abs_error_d, da_rmse_d);
        png_vis(std::string("da-") + projections[method].name, w, h, dist_area_d, true, 1.6);
        fprintf(stderr, "  isot dist (D_I) avg | max abs | rmse : %.2f | %.2f | %.2f\n", di_mean_d, di_max_abs_error_d, di_rmse_d);
        png_vis(std::string("di-") + projections[method].name, w, h, dist_isot_d, true, 1.6);

        // to copy into the LaTeX table:
#if 0
        fprintf(stderr, "  & ~\\\\~&~&~& %.3f | %.3f | %.3f & %.3f | %.3f | %.3f & ~\\\\\\hline\n",
                da_mean_d, da_max_abs_error_d, da_rmse_d,
                di_mean_d, di_max_abs_error_d, di_rmse_d);
#endif
    }
    return 0;
}
