#ifndef PROJ_OUTERRA_H
#define PROJ_OUTERRA_H

#include <type_traits>

#include "spherecube.hpp"

namespace Projection {

template <typename FLOAT>
class Outerra : public Implementation<FLOAT>
{
private:
    static constexpr FLOAT M = 1 / (2 * m_sqrt2<FLOAT>() - 2) - 1;
    int _forward_iterations;

public:
    Outerra();
    static void forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y);
    static void inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon);
};


template <typename FLOAT>
Outerra<FLOAT>::Outerra()
{
    _forward_iterations = 1;
    if (std::is_same<FLOAT, float>::value)
        _forward_iterations = 4;
    else if (std::is_same<FLOAT, double>::value)
        _forward_iterations = 5;
    else if (std::is_same<FLOAT, long double>::value)
        _forward_iterations = 5;
}

template <typename FLOAT>
void Outerra<FLOAT>::forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y)
{
    any_face_to_front_face(ctx->center, lat, lon);
    FLOAT tx = sin(lon) * cos(lat);
    FLOAT txtx = tx * tx;
    FLOAT ty = sin(lat);
    FLOAT tyty = ty * ty;
    FLOAT tz = sqrt(max(1 - txtx - tyty, static_cast<FLOAT>(0)));
    FLOAT a = M * txtx * tyty;
    FLOAT b = -M * (txtx + tyty);
    FLOAT c = -tz;
    FLOAT d = 1 + M;
    FLOAT dF;
    for (int i = 0; i < reinterpret_cast<const Outerra<FLOAT>*>(imp)->_forward_iterations; i++) {
        FLOAT tztz = tz * tz;
        FLOAT F = a * tztz * tztz + b * tztz + c * tz + d;
        FLOAT Fp = 4 * a * tztz * tz + 2 * b * tz + c;
        dF = F / Fp;
        tz = tz - dF;
    }
    x = tz * tx;
    y = tz * ty;
}

template <typename FLOAT>
void Outerra<FLOAT>::inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* /* imp */, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon)
{
    FLOAT z = 1 + M * (1 - x * x) * (1 - y * y);
    FLOAT l = sqrt(x * x + y * y + z * z);
    lat = asin(y / l);
    lon = asin((x / l) / cos(lat));
    any_face_from_front_face(ctx->center, lat, lon);
}

}

#endif
