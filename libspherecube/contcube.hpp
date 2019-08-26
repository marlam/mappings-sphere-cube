#ifndef PROJ_CONTCUBE_H
#define PROJ_CONTCUBE_H

#include "spherecube.hpp"

namespace Projection {

template <typename FLOAT>
class ContCube : public Implementation<FLOAT>
{
private:
    FLOAT C;

public:
    ContCube();
    static void forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y);
    static void inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon);
};


template <typename FLOAT>
ContCube<FLOAT>::ContCube()
{
    C = asin(1 / sqrt(static_cast<FLOAT>(3)));
}

template <typename FLOAT>
void ContCube<FLOAT>::forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y)
{
#if 0
    /* VARIANT 1: FRONT FACE, based on lat/lon */
    /* Rotate into the front cube face. Note that the original paper uses
     * a formulation for the right side cube face. */
    any_face_to_front_face(ctx->center, lat, lon);
    FLOAT cot = 1 / tan(lon + m_pi_2<FLOAT>());
    x = -atan(m_sqrt1_2<FLOAT>() * cot) / C;
    y = lat / asin(1 / sqrt(2 + cot * cot));
#endif
#if 0
    /* VARIANT 2: RIGHT FACE, based on lat/lon */
    any_face_to_front_face(ctx->center, lat, lon);
    side_face_from_front_face(CenterEquatorRight, lon);
    FLOAT cot = 1 / tan(lon);
    x = -atan(m_sqrt1_2<FLOAT>() * cot) / C;
    y = lat / asin(1 / sqrt(2 + cot * cot));
#endif
#if 1
    /* VARIANT 3: FRONT FACE, based on Gnomonic */
    FLOAT C = reinterpret_cast<const ContCube<FLOAT>*>(imp)->C;
    any_face_to_front_face(ctx->center, lat, lon);
    FLOAT u = tan(lon);
    FLOAT v = tan(lat) / cos(lon);
    x = atan(m_sqrt1_2<FLOAT>() * u) / C;
    y = atan(v*cos(atan(u))) / asin(1 / sqrt(2 + u * u));
#endif
}

template <typename FLOAT>
void ContCube<FLOAT>::inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon)
{
#if 0
    /* VARIANT 1: FRONT FACE, based on lat/lon */
    FLOAT cot = m_sqrt2<FLOAT>() * tan(-x * C);
    lon = atan2(1, cot) - m_pi_2<FLOAT>();
    lat = y * asin(1 / sqrt(2 + cot * cot));
    /* Rotate back into original cube face. */
    any_face_from_front_face(ctx->center, lat, lon);
#endif
#if 0
    /* VARIANT 2: RIGHT FACE, based on lat/lon */
    FLOAT cot = m_sqrt2<FLOAT>() * tan(-x * C);
    lon = atan2(1, cot);
    lat = y * asin(1 / sqrt(2 + cot * cot));
    side_face_to_front_face(CenterEquatorRight, lon);
    any_face_from_front_face(ctx->center, lat, lon);
#endif
#if 1
    /* VARIANT 3: FRONT FACE, based on Gnomonic */
    FLOAT C = reinterpret_cast<const ContCube<FLOAT>*>(imp)->C;
    FLOAT u = x;
    FLOAT v = y;
    x = m_sqrt2<FLOAT>() * tan(u * C);
    y = tan(v * asin(1 / sqrt(2 + x * x))) / cos(atan(x));
    lon = atan(x);
    lat = atan(y * cos(lon));
    any_face_from_front_face(ctx->center, lat, lon);
#endif
}

}

#endif
