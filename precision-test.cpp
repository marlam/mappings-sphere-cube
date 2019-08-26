#include <cmath>
#include <cstdio>
#include <type_traits>

#include "spherecube.hpp"

#include "common.hpp"


template <typename FLOAT>
void latlon_to_sphere(FLOAT lat, FLOAT lon, FLOAT& u, FLOAT& v, FLOAT& w)
{
    u = std::cos(lat) * std::cos(lon);
    v = std::cos(lat) * std::sin(lon);
    w = std::sin(lat);
}

template <typename FLOAT>
FLOAT dist_rad_sphere(FLOAT lat0, FLOAT lon0, FLOAT lat1, FLOAT lon1)
{
    /* angle between vectors: not precise enough!
    FLOAT u0, v0, w0, u1, v1, w1;
    latlon_to_sphere(lat0, lon0, u0, v0, w0);
    latlon_to_sphere(lat1, lon1, u1, v1, w1);
    FLOAT angle = std::acos(u0 * u1 + v0 * v1 + w0 * w1);
    return angle;
    */
    FLOAT slat2 = std::sin(std::fabs(lat0 - lat1) / 2);
    slat2 *= slat2;
    FLOAT slon2 = std::sin(std::fabs(lon0 - lon1) / 2);
    slon2 *= slon2;
    FLOAT t = slat2 + std::cos(lat0) * std::cos(lat1) * slon2;
    return 2.0 * std::asin(std::sqrt(t));
}

template <typename FLOAT>
FLOAT check_precision(Context<FLOAT>& ctx)
{
    const int xsteps = 2000;
    const int ysteps = 2000;

    FLOAT max_error = -1;
    int max_considered_points = 0;
    int considered_points = 0;
    for (int iy = 0; iy < ysteps; iy++) {
        FLOAT y = (iy / (ysteps - static_cast<FLOAT>(1))) * 2 - 1;
        for (int ix = 0; ix < xsteps; ix++) {
            FLOAT x = (ix / (xsteps - static_cast<FLOAT>(1))) * 2 - 1;
            FLOAT lat, lon;
            max_considered_points++;
            ctx.inverse(x, y, lat, lon);
            if (!std::isfinite(lat) || !std::isfinite(lon)) {
                //fprintf(stderr, "!cannot inverse map %g,%g! ", x, y);
                return NAN;
            }
            FLOAT fx, fy, flat, flon;
            ctx.forward(lat, lon, fx, fy);
            if (!std::isfinite(fx) || !std::isfinite(fy)) {
                // we cannot compute an error for this point
                //fprintf(stderr, "!cannot map %.19g,%.19g (comes from %g,%g)! ", degrees(lat), degrees(lon), x, y);
                continue;
            }
            ctx.inverse(fx, fy, flat, flon);
            if (!std::isfinite(flat) || !std::isfinite(flon)) {
                //fprintf(stderr, "!cannot inverse map %g,%g! ", fx, fy);
                return NAN;
            }
            considered_points++;
            FLOAT error = dist_rad_sphere(lat, lon, flat, flon);
            /*
            if (std::is_same<FLOAT, double>::value && error > 1e-3) {
                fprintf(stderr, "suspiciously large error at "
                        "x=%.17g y=%.17g lat=%.17g lon=%.17g "
                        "fx=%.17g fy=%.17g flat=%.17g flon=%.17g\n",
                        x, y, degrees(lat), degrees(lon),
                        fx, fy, degrees(flat), degrees(flon));
            }
            */
            if (error > max_error) {
                max_error = error;
            }
        }
    }
    if (considered_points < max_considered_points - max_considered_points / 100) {
        fprintf(stderr, "!only %d out of %d points could be considered! ", considered_points, max_considered_points);
        return NAN;
    }
    return max_error;
}

template <typename FLOAT>
FLOAT in_earth_millimeters(FLOAT rad_error)
{
    return rad_error * (40075.017f * 1e3f * 1e3f);
}

/*
 * Main function.
 */

int main(void)
{
    // Loop over the projection methods
    for (int method = 0; method < n_projections; method++) {
        Context<float> ctx_f(projections[method].proj,
                projections[method].params,
                projections[method].center);
        float max_error_f = check_precision(ctx_f);
        Context<double> ctx_d(projections[method].proj,
                projections[method].params,
                projections[method].center);
        double max_error_d = check_precision(ctx_d);
        /*
        Context<long double> ctx_ld(projections[method].proj,
                projections[method].params,
                projections[method].center);
        long double max_error_ld = check_precision(ctx_ld);
        */
        fprintf(stderr, "%s: %g/%g (ca. %g/%g mm on Earth)\n", projections[method].name,
                max_error_f, max_error_d,
                in_earth_millimeters(max_error_f), in_earth_millimeters(max_error_d));
        /*
        double lat, lon, x, y;
        fprintf(stderr, "%s:\n", projections[method].name);
        lat = radians(90.0);
        lon = radians(0.0);
        ctx_d.forward(lat, lon, x, y);
        fprintf(stderr, "  lat=%g lon=%g  =>   %g %g\n", degrees(lat), degrees(lon), x, y);
        lat = radians(45.0);
        lon = radians(0.0);
        ctx_d.forward(lat, lon, x, y);
        fprintf(stderr, "  lat=%g lon=%g  =>   %g %g\n", degrees(lat), degrees(lon), x, y);
        lat = radians(45.0);
        lon = radians(90.0);
        ctx_d.forward(lat, lon, x, y);
        fprintf(stderr, "  lat=%g lon=%g  =>   %g %g\n", degrees(lat), degrees(lon), x, y);
        lat = radians(45.0);
        lon = radians(180.0);
        ctx_d.forward(lat, lon, x, y);
        fprintf(stderr, "  lat=%g lon=%g  =>   %g %g\n", degrees(lat), degrees(lon), x, y);
        lat = radians(45.0);
        lon = radians(-90.0);
        ctx_d.forward(lat, lon, x, y);
        fprintf(stderr, "  lat=%g lon=%g  =>   %g %g\n", degrees(lat), degrees(lon), x, y);
        */
    }
    return 0;
}
