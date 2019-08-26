#ifndef PROJ_GNOM_TANH_H
#define PROJ_GNOM_TANH_H

#include "spherecube.hpp"

namespace Projection {

template <typename FLOAT>
class GnomonicTanh : public Implementation<FLOAT>
{
public:
    FLOAT c;
    GnomonicTanh(FLOAT c);
    static void forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y);
    static void inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon);
};


template <typename FLOAT>
GnomonicTanh<FLOAT>::GnomonicTanh(FLOAT cval) : c(cval)
{
}

template <typename FLOAT>
void GnomonicTanh<FLOAT>::forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y)
{
    any_face_to_front_face(ctx->center, lat, lon);
    x = tan(lon);
    y = tan(lat) / cos(lon);

    FLOAT c = reinterpret_cast<const GnomonicTanh<FLOAT>*>(imp)->c;
    x = 1 / tanh(c) * tanh(c * x);
    y = 1 / tanh(c) * tanh(c * y);
}

template <typename FLOAT>
void GnomonicTanh<FLOAT>::inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon)
{
    FLOAT c = reinterpret_cast<const GnomonicTanh<FLOAT>*>(imp)->c;
    x = atanh(x * tanh(c)) / c;
    y = atanh(y * tanh(c)) / c;

    lon = atan(x);
    lat = atan(y * cos(lon));
    any_face_from_front_face(ctx->center, lat, lon);
}

}

#endif
