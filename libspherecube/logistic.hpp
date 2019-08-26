#ifndef PROJ_GNOM_LOGISTIC_H
#define PROJ_GNOM_LOGISTIC_H

#include "spherecube.hpp"

namespace Projection {

template <typename FLOAT>
class GnomonicLogistic : public Implementation<FLOAT>
{
public:
    FLOAT c;
    GnomonicLogistic(FLOAT c);
    static void forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y);
    static void inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon);
};


template <typename FLOAT>
GnomonicLogistic<FLOAT>::GnomonicLogistic(FLOAT cval) : c(cval)
{
}

template <typename FLOAT>
void GnomonicLogistic<FLOAT>::forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y)
{
    any_face_to_front_face(ctx->center, lat, lon);
    x = tan(lon);
    y = tan(lat) / cos(lon);

    FLOAT c = reinterpret_cast<const GnomonicLogistic<FLOAT>*>(imp)->c;
    FLOAT d = (2*exp(c)+2)/(exp(c)-1);
    //equivalent to the above: FLOAT d = 1 / (1 / (1 + exp(-c)) - 0.5);
    x = d * (1 / (1 + exp(-c*x)) - 0.5);
    y = d * (1 / (1 + exp(-c*y)) - 0.5);
}

template <typename FLOAT>
void GnomonicLogistic<FLOAT>::inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon)
{
    FLOAT c = reinterpret_cast<const GnomonicLogistic<FLOAT>*>(imp)->c;
    x = log(-(exp(c) * x + exp(c) - x + 1) / (exp(c) * x - exp(c) - x - 1)) / c;
    y = log(-(exp(c) * y + exp(c) - y + 1) / (exp(c) * y - exp(c) - y - 1)) / c;

    lon = atan(x);
    lat = atan(y * cos(lon));
    any_face_from_front_face(ctx->center, lat, lon);
}

}

#endif
