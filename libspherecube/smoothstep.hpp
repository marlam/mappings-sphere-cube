#ifndef PROJ_GNOM_SMOOTHSTEP_H
#define PROJ_GNOM_SMOOTHSTEP_H

#include "spherecube.hpp"

namespace Projection {

template <typename FLOAT>
class GnomonicSmoothstep : public Implementation<FLOAT>
{
private:
    FLOAT c;
public:
    GnomonicSmoothstep(FLOAT cval);
    static void forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y);
    static void inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon);
};


template <typename FLOAT>
GnomonicSmoothstep<FLOAT>::GnomonicSmoothstep(FLOAT cval) : c(cval)
{
}

template <typename FLOAT>
void GnomonicSmoothstep<FLOAT>::forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y)
{
    any_face_to_front_face(ctx->center, lat, lon);
    x = tan(lon);
    y = tan(lat) / cos(lon);

    FLOAT c = reinterpret_cast<const GnomonicSmoothstep<FLOAT>*>(imp)->c;
    x = (static_cast<FLOAT>(1)/2 * (3*(c*x+1)*(c*x+1) - (c*x+1)*(c*x+1)*(c*x+1)) -1)
        / (static_cast<FLOAT>(1)/2 * (3*(c+1)*(c+1) - (c+1)*(c+1)*(c+1)) -1);
    y = (static_cast<FLOAT>(1)/2 * (3*(c*y+1)*(c*y+1) - (c*y+1)*(c*y+1)*(c*y+1)) -1)
        / (static_cast<FLOAT>(1)/2 * (3*(c+1)*(c+1) - (c+1)*(c+1)*(c+1)) -1);
}

template <typename FLOAT>
void GnomonicSmoothstep<FLOAT>::inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon)
{
    FLOAT c = reinterpret_cast<const GnomonicSmoothstep<FLOAT>*>(imp)->c;
    x = (static_cast<FLOAT>(-1)/2 * (2-c*x)*(c*x+1)*(c*x+1) + 2*c*x + 1)
        / (static_cast<FLOAT>(-1)/2 * (2-c)*(c+1)*(c+1) + 2*c + 1);
    y = (static_cast<FLOAT>(-1)/2 * (2-c*y)*(c*y+1)*(c*y+1) + 2*c*y + 1)
        / (static_cast<FLOAT>(-1)/2 * (2-c)*(c+1)*(c+1) + 2*c + 1);

    lon = atan(x);
    lat = atan(y * cos(lon));
    any_face_from_front_face(ctx->center, lat, lon);
}

}

#endif
