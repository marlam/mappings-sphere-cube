#ifndef PROJ_GNOM_EVERITT_H
#define PROJ_GNOM_EVERITT_H

#include "spherecube.hpp"

namespace Projection {

template <typename FLOAT>
class GnomonicEveritt : public Implementation<FLOAT>
{
public:
    FLOAT c;
    GnomonicEveritt(FLOAT c);
    static void forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y);
    static void inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon);
};


template <typename FLOAT>
GnomonicEveritt<FLOAT>::GnomonicEveritt(FLOAT cval) : c(cval)
{
}

template <typename FLOAT>
void GnomonicEveritt<FLOAT>::forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y)
{
    any_face_to_front_face(ctx->center, lat, lon);
    x = tan(lon);
    y = tan(lat) / cos(lon);

    FLOAT c = reinterpret_cast<const GnomonicEveritt<FLOAT>*>(imp)->c;
    x = x * (c + (1 - c) * fabs(x));
    y = y * (c + (1 - c) * fabs(y));
}

template <typename FLOAT>
void GnomonicEveritt<FLOAT>::inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon)
{
    FLOAT c = reinterpret_cast<const GnomonicEveritt<FLOAT>*>(imp)->c;
    x = sign_not_zero(x) * (c - sqrt(c * c - 4 * (c - 1) * fabs(x))) / (2 * (c - 1));
    y = sign_not_zero(y) * (c - sqrt(c * c - 4 * (c - 1) * fabs(y))) / (2 * (c - 1));

    lon = atan(x);
    lat = atan(y * cos(lon));
    any_face_from_front_face(ctx->center, lat, lon);
}

}

#endif
