#ifndef PROJ_GNOM_ATAN_H
#define PROJ_GNOM_ATAN_H

#include "spherecube.hpp"

namespace Projection {

template <typename FLOAT>
class GnomonicTangens : public Implementation<FLOAT>
{
public:
    FLOAT c;
    GnomonicTangens(FLOAT c);
    static void forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y);
    static void inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon);
};


template <typename FLOAT>
GnomonicTangens<FLOAT>::GnomonicTangens(FLOAT cval) : c(cval)
{
}

template <typename FLOAT>
void GnomonicTangens<FLOAT>::forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y)
{
    any_face_to_front_face(ctx->center, lat, lon);
    x = tan(lon);
    y = tan(lat) / cos(lon);

    FLOAT c = reinterpret_cast<const GnomonicTangens<FLOAT>*>(imp)->c;
    x = 1 / atan(c) * atan(c * x);
    y = 1 / atan(c) * atan(c * y);

}

template <typename FLOAT>
void GnomonicTangens<FLOAT>::inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon)
{
    FLOAT c = reinterpret_cast<const GnomonicTangens<FLOAT>*>(imp)->c;
    x = tan(x * atan(c)) / c;
    y = tan(y * atan(c)) / c;

    lon = atan(x);
    lat = atan(y * cos(lon));
    any_face_from_front_face(ctx->center, lat, lon);
}

}

#endif
