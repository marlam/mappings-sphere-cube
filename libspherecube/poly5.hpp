#ifndef PROJ_GNOM_POLY5_H
#define PROJ_GNOM_POLY5_H

#include "spherecube.hpp"

namespace Projection {

template <typename FLOAT>
class GnomonicPoly5 : public Implementation<FLOAT>
{
public:
    const bool iterativeInverse;
    GnomonicPoly5(bool ii);
    static void forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y);
    static void inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon);
};


/* From: "Cube-to-sphere Projections for Procedural Texturing and Beyond",
 * Zucker and Higashi, JCGT vol 7 no 2, 2018
 *
 * But note that the use of parameters differ from their description!
 * This can be seen in their source code. */

template <typename FLOAT>
GnomonicPoly5<FLOAT>::GnomonicPoly5(bool ii) : iterativeInverse(ii)
{
}

template <typename FLOAT> FLOAT func(FLOAT a)
{
    FLOAT a2 = a * a;
    return (0.745558715593 + (0.130546850193 + 0.123894434214 * a2) * a2) * a;
}

template <typename FLOAT> FLOAT inv_func(FLOAT a)
{
    FLOAT a2 = a * a;
    return (1.34318229552 + (-0.486514066449 + 0.143331770927 * a2) * a2) * a;
}

template <typename FLOAT>
void GnomonicPoly5<FLOAT>::forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y)
{
    any_face_to_front_face(ctx->center, lat, lon);
    x = tan(lon);
    y = tan(lat) / cos(lon);

    FLOAT a = inv_func(x);
    FLOAT b = inv_func(y);

    if (reinterpret_cast<const GnomonicPoly5<FLOAT>*>(imp)->iterativeInverse) {
        FLOAT err0_a = func(a) - x;
        FLOAT err0_b = func(b) - y;
        FLOAT delta_a = static_cast<FLOAT>(-0.1) * err0_a;  // -0.1 is used by Zucker; -1 actually works better here
        FLOAT delta_b = static_cast<FLOAT>(-0.1) * err0_b;  // -0.1 is used by Zucker; -1 actually works better here
        a += delta_a;
        b += delta_b;
        for (int i = 0; i < 9; i++) {
            FLOAT err1_a = func(a) - x;
            FLOAT err1_b = func(b) - y;
            FLOAT delta_err_a = err1_a - err0_a;
            FLOAT delta_err_b = err1_b - err0_a;
            if (abs(delta_err_a) > static_cast<FLOAT>(1e-15)) { // unclear if 1e-15 is ok; Zucker uses 1e-7 for float
                FLOAT inv_slope_a = delta_a / delta_err_a;
                delta_a = -err1_a * inv_slope_a;
                a += delta_a;
                err0_a = err1_a;
            }
            if (abs(delta_err_b) > static_cast<FLOAT>(1e-15)) { // unclear if 1e-15 is ok; Zucker uses 1e-7 for float
                FLOAT inv_slope_b = delta_b / delta_err_b;
                delta_b = -err1_b * inv_slope_b;
                b += delta_b;
                err0_b = err1_b;
            }
        }
    }

    x = a;
    y = b;
}

template <typename FLOAT>
void GnomonicPoly5<FLOAT>::inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>*, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon)
{
    x = func(x);
    y = func(y);

    lon = atan(x);
    lat = atan(y * cos(lon));
    any_face_from_front_face(ctx->center, lat, lon);
}

}

#endif
