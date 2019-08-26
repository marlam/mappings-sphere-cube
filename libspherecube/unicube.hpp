#ifndef PROJ_UNICUBE_H
#define PROJ_UNICUBE_H

#include "spherecube.hpp"

namespace Projection {

template <typename FLOAT>
class UniCube : public Implementation<FLOAT>
{
public:
    UniCube();
    static void forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y);
    static void inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon);
};


template <typename FLOAT>
UniCube<FLOAT>::UniCube()
{
}

template <typename FLOAT>
void UniCube<FLOAT>::forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* /* imp */, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y)
{
    // Rotate into the front cube face
    any_face_to_front_face(ctx->center, lat, lon);
    // Gnomonic projection
    x = tan(lon);
    y = tan(lat) / cos(lon);
    // UniCube adjustment, see Eq. (1) in paper
    x = 6 / m_pi<FLOAT>() * asin(x / sqrt(2 * x * x + 2));
    y = 6 / m_pi<FLOAT>() * asin(y / sqrt(2 * y * y + 2));
}

template <typename FLOAT>
void UniCube<FLOAT>::inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* /* imp */, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon)
{
    // Undo the UniCube adjustment, see Eq. (13) in paper
    FLOAT sx = sin(m_pi<FLOAT>() / 6 * x);
    x = sx / sqrt(static_cast<FLOAT>(0.5) - sx * sx);
    FLOAT sy = sin(m_pi<FLOAT>() / 6 * y);
    y = sy / sqrt(static_cast<FLOAT>(0.5) - sy * sy);
    // Inverse Gnomonic projection
    lon = atan(x);
    lat = atan(y * cos(lon));
    // Rotate back into original cube face
    any_face_from_front_face(ctx->center, lat, lon);
}

}

#endif
