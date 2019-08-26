#ifndef PROJ_NOWELL_H
#define PROJ_NOWELL_H

#include "spherecube.hpp"

namespace Projection {

template <typename FLOAT>
class Nowell : public Implementation<FLOAT>
{
public:
    Nowell();
    static void forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y);
    static void inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon);
};

/*
 * This is based on a mapping between unit sphere and unit cube,
 * developed and described in various internet sources, most importantly
 * the following:
 *
 * - Original idea, sphere -> cube transformation:
 *   http://mathproofs.blogspot.de/2005/07/mapping-cube-to-sphere.html
 *
 * - Inverse transformation cube -> sphere:
 *   http://stackoverflow.com/questions/2656899/mapping-a-sphere-to-a-cube
 */

template <typename FLOAT>
Nowell<FLOAT>::Nowell()
{
}

template <typename FLOAT>
void Nowell<FLOAT>::forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* /* imp */, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y)
{
    FLOAT sinlat = sin(lat);
    FLOAT coslat = cos(lat);
    FLOAT sinlon = sin(lon);
    FLOAT coslon = cos(lon);

    // This is to avoid a precision problem with lon=pi.
    // sin(pi) and cos(pi) are not what they should be,
    // and this error propagates in an unfortunate manner
    // so that the resulting error of forward and inverse
    // projection on an Earth-sized planet can reach half
    // a meter. With this workaround, we are far in the
    // sub-millimeter range.
    if (fabs(lon) >= m_pi<FLOAT>()) {
        sinlon = 0;
        coslon = -1;
    }

    // Get cartesian coordinates on the sphere.
    FLOAT u = coslat * coslon;
    FLOAT v = coslat * sinlon;
    FLOAT w = sinlat;

    // Transform to cube coordinates, using code adapted from here:
    // http://stackoverflow.com/questions/2656899/mapping-a-sphere-to-a-cube
    if (ctx->center == CenterEquatorLeft || ctx->center == CenterEquatorRight) {
        FLOAT aa2 = u * u * 2;
        FLOAT bb2 = w * w * 2;
        FLOAT inner = bb2 - aa2 - 3;
        FLOAT innersqrt = -sqrt(inner * inner - 12 * aa2);
        x = sign(u) * sqrt(max(innersqrt + aa2 - bb2 + 3, static_cast<FLOAT>(0))) * m_sqrt1_2<FLOAT>();
        y = sign(w) * sqrt(max(innersqrt - aa2 + bb2 + 3, static_cast<FLOAT>(0))) * m_sqrt1_2<FLOAT>();
        if (ctx->center == CenterEquatorRight)
            x = -x;
    } else if (ctx->center == CenterEquatorFront || ctx->center == CenterEquatorBack) {
        FLOAT aa2 = v * v * 2;
        FLOAT bb2 = w * w * 2;
        FLOAT inner = bb2 - aa2 - 3;
        FLOAT innersqrt = -sqrt(inner * inner - 12 * aa2);
        x = sign(v) * sqrt(max(innersqrt + aa2 - bb2 + 3, static_cast<FLOAT>(0))) * m_sqrt1_2<FLOAT>();
        y = sign(w) * sqrt(max(innersqrt - aa2 + bb2 + 3, static_cast<FLOAT>(0))) * m_sqrt1_2<FLOAT>();
        if (ctx->center == CenterEquatorBack)
            x = -x;
    } else /* CenterNorthPole or CenterSouthPole */ {
        FLOAT aa2 = u * u * 2;
        FLOAT bb2 = v * v * 2;
        FLOAT inner = bb2 - aa2 - 3;
        FLOAT innersqrt = -sqrt(inner * inner - 12 * aa2);
        y = sign(u) * sqrt(max(innersqrt + aa2 - bb2 + 3, static_cast<FLOAT>(0))) * m_sqrt1_2<FLOAT>();
        x = sign(v) * sqrt(max(innersqrt - aa2 + bb2 + 3, static_cast<FLOAT>(0))) * m_sqrt1_2<FLOAT>();
        if (ctx->center == CenterNorthPole)
            y = -y;
    }
}

template <typename FLOAT>
void Nowell<FLOAT>::inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* /* imp */, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon)
{
    FLOAT cx, cy, cz;
    FLOAT xx, yy, zz;
    FLOAT u, v, w;

    // Get cube coordinates from map coordinates. Depends on the cube side.
    if (ctx->center == CenterEquatorFront) {
        cx = +1;
        cy = x;
        cz = y;
    } else if (ctx->center == CenterEquatorRight) {
        cx = -x;
        cy = +1;
        cz = y;
    } else if (ctx->center == CenterEquatorBack) {
        cx = -1;
        cy = -x;
        cz = y;
    } else if (ctx->center == CenterEquatorLeft) {
        cx = x;
        cy = -1;
        cz = y;
    } else if (ctx->center == CenterNorthPole) {
        cx = -y;
        cy = x;
        cz = +1;
    } else /* CenterSouthPole */ {
        cx = y;
        cy = x;
        cz = -1;
    }

    // Compute sphere coordinates according to
    // http://mathproofs.blogspot.de/2005/07/mapping-cube-to-sphere.html
    xx = cx * cx;
    yy = cy * cy;
    zz = cz * cz;
    u = cx * sqrt(1 - yy / 2 - zz / 2 + yy * zz / 3);
    v = cy * sqrt(1 - xx / 2 - zz / 2 + xx * zz / 3);
    w = cz * sqrt(1 - xx / 2 - yy / 2 + xx * yy / 3);

    // Compute lat, lon from sphere coordinates
    lat = asin(w);
    lon = atan2(v, u);
}

}

#endif
