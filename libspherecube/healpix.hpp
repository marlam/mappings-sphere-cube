#ifndef PROJ_HEALPIX_H
#define PROJ_HEALPIX_H

#include "spherecube.hpp"

namespace Projection {

template <typename FLOAT>
class HEALPix : public Implementation<FLOAT>
{
public:
    HEALPix();
    static void forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y);
    static void inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon);
};


template <typename FLOAT>
HEALPix<FLOAT>::HEALPix()
{
}

template <typename FLOAT>
void HEALPix<FLOAT>::forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* /* imp */, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y)
{
    if (ctx->center == CenterNorthPole || ctx->center == CenterSouthPole) {
        int area = (fabs(lon) <= m_pi_4<FLOAT>() ? 0 : fabs(lon) >= 3 * m_pi_4<FLOAT>() ? 2 : lon > 0 ? 1 : 3);
        FLOAT s = sqrt(3 * max(1 - fabs(sin(lat)), static_cast<FLOAT>(0)));
        FLOAT c = (area > 2 ? 1 : area) * sign(lon) * m_pi_2<FLOAT>();
        FLOAT t = (lon - c) * (4 / m_pi<FLOAT>());
        FLOAT x_sign = (area == 2 || area == 3 ? -1 : +1);
        FLOAT y_sign = (area == 0 || area == 3 ? -1 : +1) * (ctx->center == CenterSouthPole ? -1 : +1);
        FLOAT x_factor = (area == 0 || area == 2 ? t : 1);
        FLOAT y_factor = (area == 0 || area == 2 ? 1 : t);
        x = x_sign * s * x_factor;
        y = y_sign * s * y_factor;
    } else {
        side_face_to_front_face(ctx->center, lon);
        x = lon * 4 / m_pi<FLOAT>();
        y = sin(lat) * 3 / 2;
    }
}

template <typename FLOAT>
void HEALPix<FLOAT>::inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* /* imp */, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon)
{
    if (ctx->center == CenterNorthPole || ctx->center == CenterSouthPole) {
        FLOAT pole_sign = (ctx->center == CenterSouthPole ? -1 : +1);
        int area = (fabs(x) <= -pole_sign * y ? 0 : fabs(x) <= pole_sign * y ? 2 : x > 0 ? 1 : 3);
        FLOAT x_sign = (area == 2 || area == 3 ? -1 : +1);
        FLOAT y_sign = (area == 0 || area == 3 ? -1 : +1) * (ctx->center == CenterSouthPole ? -1 : +1);
        FLOAT s = ((area == 0 || area == 2) ? y_sign * y : x_sign * x);
        FLOAT t = ((area == 0 || area == 2) ? x_sign * x : y_sign * y) / s;
        FLOAT c = (area > 2 ? 1 : area) * sign(x) * m_pi_2<FLOAT>();
        lat = pole_sign * asin(1 - s * s / 3);
        lon = t * m_pi_4<FLOAT>() + c;
    } else {
        lon = x * m_pi<FLOAT>() / 4;
        lat = asin(2 * y / 3);
        side_face_from_front_face(ctx->center, lon);
    }
}

}

#endif
