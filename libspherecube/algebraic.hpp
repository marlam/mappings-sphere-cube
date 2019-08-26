#ifndef PROJ_GNOM_ALGEBRAIC_H
#define PROJ_GNOM_ALGEBRAIC_H

#include "spherecube.hpp"

namespace Projection {

template <typename FLOAT>
class GnomonicAlgebraic : public Implementation<FLOAT>
{
public:
    FLOAT c;
    GnomonicAlgebraic(FLOAT c);
    static void forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y);
    static void inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon);
};


template <typename FLOAT>
GnomonicAlgebraic<FLOAT>::GnomonicAlgebraic(FLOAT cval) : c(cval)
{
}

template <typename FLOAT>
void GnomonicAlgebraic<FLOAT>::forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y)
{
    any_face_to_front_face(ctx->center, lat, lon);
    x = tan(lon);
    y = tan(lat) / cos(lon);

    FLOAT c = reinterpret_cast<const GnomonicAlgebraic<FLOAT>*>(imp)->c;
    x = x * sqrt((1+c*c) / (1 + c * x * c * x));
    y = y * sqrt((1+c*c) / (1 + c * y * c * y));
}

template <typename FLOAT>
void GnomonicAlgebraic<FLOAT>::inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon)
{
    FLOAT c = reinterpret_cast<const GnomonicAlgebraic<FLOAT>*>(imp)->c;
    x = x / sqrt(1+c*c - c*c*x*x);
    y = y / sqrt(1+c*c - c*c*y*y);

    lon = atan(x);
    lat = atan(y * cos(lon));
    any_face_from_front_face(ctx->center, lat, lon);
}

}

#endif
