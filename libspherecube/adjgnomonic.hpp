#ifndef PROJ_ADJGNOMONIC_H
#define PROJ_ADJGNOMONIC_H

#include "spherecube.hpp"

namespace Projection {

template <typename FLOAT>
class AdjGnomonic : public Implementation<FLOAT>
{
public:
    AdjGnomonic();
    static void forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y);
    static void inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon);
};


template <typename FLOAT>
AdjGnomonic<FLOAT>::AdjGnomonic()
{
}

template <typename FLOAT>
void AdjGnomonic<FLOAT>::forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* /* imp */, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y)
{
    any_face_to_front_face(ctx->center, lat, lon);
    x = 4 / m_pi<FLOAT>() * lon;
    y = 4 / m_pi<FLOAT>() * atan(tan(lat) / cos(lon));
}

template <typename FLOAT>
void AdjGnomonic<FLOAT>::inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* /* imp */, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon)
{
    lon = m_pi<FLOAT>() / 4 * x;
    lat = atan(tan(m_pi<FLOAT>() / 4 * y) * cos(lon));
    any_face_from_front_face(ctx->center, lat, lon);
}

}

#endif
