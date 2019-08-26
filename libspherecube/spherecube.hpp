#ifndef SPHERECUBE_H
#define SPHERECUBE_H

#include <cmath>
#include <string>

namespace Projection {

enum Projection {
    ProjGnomonic,
    ProjAdjGnomonic, // equivalent to ProjGnomTangens for c=1
    ProjQSC,
    ProjNowell,
    ProjOuterra,
    ProjHEALPix,
    ProjContCube,
    ProjUniCube,
    ProjGnomTangens,
    ProjGnomTanh,
    ProjGnomLogistic,
    ProjGnomAlgebraic,
    ProjGnomPoly5,
    ProjGnomSmoothstep,
    ProjGnomEveritt,
};

enum Center {
    CenterNorthPole = 0,
    CenterSouthPole = 1,
    CenterEquatorFront = 2,
    CenterEquatorRight = 3,
    CenterEquatorBack = 4,
    CenterEquatorLeft = 5
};

using std::sin;
using std::asin;
using std::cos;
using std::acos;
using std::tan;
using std::atan;
using std::atan2;
using std::fabs;
using std::hypot;
using std::sqrt;
using std::cbrt;
using std::isfinite;
using std::min;
using std::max;
using std::floor;
using std::ceil;
using std::log;
using std::exp;
using std::tanh;
using std::atanh;
using std::abs;

template <typename FLOAT>
constexpr FLOAT m_pi() { return M_PIl; }

template <typename FLOAT>
constexpr FLOAT m_pi_2() { return M_PI_2l; }

template <typename FLOAT>
constexpr FLOAT m_pi_4() { return M_PI_4l; }

template <typename FLOAT>
constexpr FLOAT m_sqrt2() { return M_SQRT2l; }

template <typename FLOAT>
constexpr FLOAT m_sqrt1_2() { return M_SQRT1_2l; }

template <typename FLOAT>
constexpr FLOAT sign(FLOAT x) { return (x < 0 ? -1 : x > 0 ? +1 : 0); }

template <typename FLOAT>
constexpr FLOAT sign_not_zero(FLOAT x) { return (x < 0 ? -1 : +1); }

template <typename FLOAT>
constexpr FLOAT degrees(FLOAT x) { return 180 / m_pi<FLOAT>() * x; }

template <typename FLOAT>
constexpr FLOAT radians(FLOAT x) { return m_pi<FLOAT>() / 180 * x; }

template <typename FLOAT>
constexpr FLOAT clamp(FLOAT x, FLOAT lo, FLOAT hi) { return x < lo ? lo : x > hi ? hi : x; }

template <typename FLOAT>
bool possibly_in_cubeside(Center c, FLOAT lat, FLOAT lon, FLOAT tolerance = radians<FLOAT>(1))
{
    // asin(1/sqrt(3)) = 0.61547970867
    return (
            c == CenterNorthPole ? (lat >= static_cast<FLOAT>(0.61547970) - tolerance)
            : c == CenterSouthPole ? (lat <= static_cast<FLOAT>(-0.61547970) + tolerance)
            : c == CenterEquatorFront ? (fabs(lat) <= m_pi_4<FLOAT>() + tolerance && fabs(lon) <= m_pi_4<FLOAT>() + tolerance)
            : c == CenterEquatorRight ? (fabs(lat) <= m_pi_4<FLOAT>() + tolerance && lon >= m_pi_4<FLOAT>() - tolerance && lon <= 3 * m_pi_4<FLOAT>() + tolerance)
            : c == CenterEquatorBack ? (fabs(lat) <= m_pi_4<FLOAT>() + tolerance && fabs(lon) >= 3 * m_pi_4<FLOAT>() - tolerance)
            : (fabs(lat) <= m_pi_4<FLOAT>() + tolerance && lon <= -m_pi_4<FLOAT>() + tolerance && lon >= -3 * m_pi_4<FLOAT>() - tolerance));
}

template <typename FLOAT>
bool in_cubeside(Center c, FLOAT sphere_x, FLOAT sphere_y, FLOAT sphere_z)
{
    return (
            c == CenterNorthPole ?      (sphere_z >= sphere_x && sphere_z >= sphere_y)
            : c == CenterSouthPole ?    (sphere_z <= sphere_x && sphere_z <= sphere_y)
            : c == CenterEquatorFront ? (sphere_x >= sphere_y && sphere_x >= sphere_z)
            : c == CenterEquatorRight ? (sphere_y >= sphere_x && sphere_y >= sphere_z)
            : c == CenterEquatorBack ?  (sphere_x <= sphere_y && sphere_x <= sphere_z)
            :                           (sphere_y <= sphere_x && sphere_y <= sphere_z));
}

template <typename FLOAT>
bool in_hemisphere(Center c, FLOAT lat, FLOAT lon)
{
    return (
            c == CenterNorthPole ?      (lat >= 0)
            : c == CenterSouthPole ?    (lat <= 0)
            : c == CenterEquatorFront ? (fabs(lon) <= m_pi_2<FLOAT>())
            : c == CenterEquatorRight ? (lon >= 0 && lon <= m_pi<FLOAT>())
            : c == CenterEquatorBack ?  (fabs(lon) >= m_pi_2<FLOAT>())
            :                           (lon <= 0 && lon >= -m_pi<FLOAT>()));
}

// rotate 90 degrees around x axis
template <typename FLOAT>
void rot_x_fwd(FLOAT& y, FLOAT& z)
{
    FLOAT tmp = y;
    y = -z;
    z = tmp;
}

// rotate -90 degrees around x axis
template <typename FLOAT>
void rot_x_bwd(FLOAT& y, FLOAT& z)
{
    FLOAT tmp = y;
    y = z;
    z = -tmp;
}

// rotate 90 degrees around y axis
template <typename FLOAT>
void rot_y_fwd(FLOAT& x, FLOAT& z)
{
    FLOAT tmp = x;
    x = z;
    z = -tmp;
}

// rotate -90 degrees around y axis
template <typename FLOAT>
void rot_y_bwd(FLOAT& x, FLOAT& z)
{
    FLOAT tmp = x;
    x = -z;
    z = tmp;
}

// rotate 90 degrees around z axis
template <typename FLOAT>
void rot_z_fwd(FLOAT& x, FLOAT& y)
{
    FLOAT tmp = x;
    x = -y;
    y = tmp;
}

// rotate -90 degrees around z axis
template <typename FLOAT>
void rot_z_bwd(FLOAT& x, FLOAT& y)
{
    FLOAT tmp = x;
    x = y;
    y = -tmp;
}

// convert geodetic lat, lon to cartesian x, y, z
template <typename FLOAT>
void to_cartesian(const FLOAT lat, const FLOAT lon, FLOAT& x, FLOAT& y, FLOAT& z)
{
    FLOAT coslat = cos(lat);
    FLOAT sinlat = sin(lat);
    FLOAT coslon = cos(lon);
    FLOAT sinlon = sin(lon);
    x = coslat * coslon;
    y = coslat * sinlon;
    z = sinlat;
}

// convert cartesian x, y, z to geodetic lat, lon
template <typename FLOAT>
void to_geodetic(const FLOAT x, const FLOAT y, const FLOAT z, FLOAT& lat, FLOAT& lon)
{
    lat = asin(z);
    lon = atan2(y, x);
}

// rotate the given side face (front/back/left/right) to the front face
template <typename FLOAT>
void side_face_to_front_face(Center c, FLOAT& lon)
{
    if (c == CenterEquatorRight) {
        lon -= m_pi_2<FLOAT>();
        if (lon < -m_pi<FLOAT>())
            lon += 2 * m_pi<FLOAT>();
    } else if (c == CenterEquatorBack) {
        lon -= m_pi<FLOAT>();
        if (lon < -m_pi<FLOAT>())
            lon += 2 * m_pi<FLOAT>();
    } else if (c == CenterEquatorLeft) {
        lon += m_pi_2<FLOAT>();
        if (lon > m_pi<FLOAT>())
            lon -= 2 * m_pi<FLOAT>();
    }
}

// rotate the front face to the given side face (front/back/left/right)
template <typename FLOAT>
void side_face_from_front_face(Center c, FLOAT& lon)
{
    if (c == CenterEquatorRight) {
        lon += m_pi_2<FLOAT>();
        if (lon > m_pi<FLOAT>())
            lon -= 2 * m_pi<FLOAT>();
    } else if (c == CenterEquatorBack) {
        lon += m_pi<FLOAT>();
        if (lon > m_pi<FLOAT>())
            lon -= 2 * m_pi<FLOAT>();
    } else if (c == CenterEquatorLeft) {
        lon -= m_pi_2<FLOAT>();
        if (lon < -m_pi<FLOAT>())
            lon += 2 * m_pi<FLOAT>();
    }
}

// rotate any face to the front face
template <typename FLOAT>
void any_face_to_front_face(Center c, FLOAT& lat, FLOAT& lon)
{
    FLOAT x = 0, y = 0, z = 0;
    if (c == CenterNorthPole) {
        to_cartesian(lat, lon, x, y, z);
        rot_y_fwd(x, z);
        to_geodetic(x, y, z, lat, lon);
    } else if (c == CenterSouthPole) {
        to_cartesian(lat, lon, x, y, z);
        rot_y_bwd(x, z);
        to_geodetic(x, y, z, lat, lon);
    } else {
        side_face_to_front_face(c, lon);
    }
}

// rotate the front face to any face
template <typename FLOAT>
void any_face_from_front_face(Center c, FLOAT& lat, FLOAT& lon)
{
    FLOAT x = 0, y = 0, z = 0;
    if (c == CenterNorthPole) {
        to_cartesian(lat, lon, x, y, z);
        rot_y_bwd(x, z);
        to_geodetic(x, y, z, lat, lon);
    } else if (c == CenterSouthPole) {
        to_cartesian(lat, lon, x, y, z);
        rot_y_fwd(x, z);
        to_geodetic(x, y, z, lat, lon);
    } else {
        side_face_from_front_face(c, lon);
    }
}

// rotate any face to the top face
template <typename FLOAT>
void any_face_to_top_face(Center c, FLOAT& lat, FLOAT& lon)
{
    if (c == CenterNorthPole) {
        // do nothing
    } else if (c == CenterSouthPole) {
        lat = -lat;
        lon = sign_not_zero(lon) * m_pi<FLOAT>() - lon;
    } else {
        FLOAT x = 0, y = 0, z = 0;
        if (c == CenterEquatorFront) {
            to_cartesian(lat, lon, x, y, z);
            rot_y_bwd(x, z);
            to_geodetic(x, y, z, lat, lon);
        } else if (c == CenterEquatorRight) {
            to_cartesian(lat, lon, x, y, z);
            rot_x_fwd(y, z);
            to_geodetic(x, y, z, lat, lon);
            lon -= m_pi_2<FLOAT>();
            if (lon < -m_pi<FLOAT>())
                lon += 2 * m_pi<FLOAT>();
        } else if (c == CenterEquatorBack) {
            to_cartesian(lat, lon, x, y, z);
            rot_y_fwd(x, z);
            to_geodetic(x, y, z, lat, lon);
            lon += m_pi<FLOAT>();
            if (lon > m_pi<FLOAT>())
                lon -= 2 * m_pi<FLOAT>();
        } else {
            to_cartesian(lat, lon, x, y, z);
            rot_x_bwd(y, z);
            to_geodetic(x, y, z, lat, lon);
            lon += m_pi_2<FLOAT>();
            if (lon > m_pi<FLOAT>())
                lon -= 2 * m_pi<FLOAT>();
        }
    }
}

// rotate the top face to any face
template <typename FLOAT>
void any_face_from_top_face(Center c, FLOAT& lat, FLOAT& lon)
{
    if (c == CenterNorthPole) {
        // do nothing
    } else if (c == CenterSouthPole) {
        lat = -lat;
        lon = sign_not_zero(lon) * m_pi<FLOAT>() - lon;
    } else {
        FLOAT x = 0, y = 0, z = 0;
        if (c == CenterEquatorFront) {
            to_cartesian(lat, lon, x, y, z);
            rot_y_fwd(x, z);
            to_geodetic(x, y, z, lat, lon);
        } else if (c == CenterEquatorRight) {
            lon += m_pi_2<FLOAT>();
            if (lon > m_pi<FLOAT>())
                lon -= 2 * m_pi<FLOAT>();
            to_cartesian(lat, lon, x, y, z);
            rot_x_bwd(y, z);
            to_geodetic(x, y, z, lat, lon);
        } else if (c == CenterEquatorBack) {
            lon -= m_pi<FLOAT>();
            if (lon < -m_pi<FLOAT>())
                lon += 2 * m_pi<FLOAT>();
            to_cartesian(lat, lon, x, y, z);
            rot_y_bwd(x, z);
            to_geodetic(x, y, z, lat, lon);
        } else {
            lon -= m_pi_2<FLOAT>();
            if (lon < -m_pi<FLOAT>())
                lon += 2 * m_pi<FLOAT>();
            to_cartesian(lat, lon, x, y, z);
            rot_x_fwd(y, z);
            to_geodetic(x, y, z, lat, lon);
        }
    }
}

template <typename FLOAT>
class Context;

template <typename FLOAT>
class Implementation
{
public:
    static void forward(const Context<FLOAT>* /* ctx */, const Implementation<FLOAT>* /* imp */, FLOAT /* lat */, FLOAT /* lon */, FLOAT& /* x */, FLOAT& /* y */) {}
    static void inverse(const Context<FLOAT>* /* ctx */, const Implementation<FLOAT>* /* imp */, FLOAT /* x */, FLOAT /* y */, FLOAT& /* lat */, FLOAT& /* lon */) {}
};

}

#include "gnomonic.hpp"
#include "adjgnomonic.hpp"
#include "qsc.hpp"
#include "nowell.hpp"
#include "outerra.hpp"
#include "healpix.hpp"
#include "contcube.hpp"
#include "unicube.hpp"
#include "tangens.hpp"
#include "tanh.hpp"
#include "logistic.hpp"
#include "algebraic.hpp"
#include "poly5.hpp"
#include "smoothstep.hpp"
#include "everitt.hpp"

namespace Projection {

template <typename FLOAT>
class Context
{
private:
    bool _valid;
    Implementation<FLOAT>* _imp;
    void (*_forward)(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y);
    void (*_inverse)(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon);

public:

    /* Generic parameters of the projection */
    const enum Projection projection;
    const enum Center center;

    /* User interface */
    Context(enum Projection P,
            const std::string& params,
            enum Center center) :
        _valid(false),
        _imp(nullptr), _forward(nullptr), _inverse(nullptr),
        projection(P), center(center)
    {
        float param;
        switch (P) {
        case ProjGnomonic:
            _imp = new Gnomonic<FLOAT>();
            _forward = Gnomonic<FLOAT>::forward;
            _inverse = Gnomonic<FLOAT>::inverse;
            _valid = true;
            break;
        case ProjAdjGnomonic:
            _imp = new AdjGnomonic<FLOAT>();
            _forward = AdjGnomonic<FLOAT>::forward;
            _inverse = AdjGnomonic<FLOAT>::inverse;
            _valid = true;
            break;
        case ProjQSC:
            _imp = new QSC<FLOAT>();
            _forward = QSC<FLOAT>::forward;
            _inverse = QSC<FLOAT>::inverse;
            _valid = true;
            break;
        case ProjNowell:
            _imp = new Nowell<FLOAT>();
            _forward = Nowell<FLOAT>::forward;
            _inverse = Nowell<FLOAT>::inverse;
            _valid = true;
            break;
        case ProjOuterra:
            _imp = new Outerra<FLOAT>();
            _forward = Outerra<FLOAT>::forward;
            _inverse = Outerra<FLOAT>::inverse;
            _valid = true;
            break;
        case ProjHEALPix:
            _imp = new HEALPix<FLOAT>();
            _forward = HEALPix<FLOAT>::forward;
            _inverse = HEALPix<FLOAT>::inverse;
            _valid = true;
            break;
        case ProjContCube:
            _imp = new ContCube<FLOAT>();
            _forward = ContCube<FLOAT>::forward;
            _inverse = ContCube<FLOAT>::inverse;
            _valid = true;
            break;
        case ProjUniCube:
            _imp = new UniCube<FLOAT>();
            _forward = UniCube<FLOAT>::forward;
            _inverse = UniCube<FLOAT>::inverse;
            _valid = true;
            break;
        case ProjGnomTangens:
            param = 1.0f;
            sscanf(params.c_str(), "c=%f", &param);
            _imp = new GnomonicTangens<FLOAT>(param);
            _forward = GnomonicTangens<FLOAT>::forward;
            _inverse = GnomonicTangens<FLOAT>::inverse;
            _valid = true;
            break;
        case ProjGnomTanh:
            param = 1.0f;
            sscanf(params.c_str(), "c=%f", &param);
            _imp = new GnomonicTanh<FLOAT>(param);
            _forward = GnomonicTanh<FLOAT>::forward;
            _inverse = GnomonicTanh<FLOAT>::inverse;
            _valid = true;
            break;
        case ProjGnomLogistic:
            param = 1.0f;
            sscanf(params.c_str(), "c=%f", &param);
            _imp = new GnomonicLogistic<FLOAT>(param);
            _forward = GnomonicLogistic<FLOAT>::forward;
            _inverse = GnomonicLogistic<FLOAT>::inverse;
            _valid = true;
            break;
        case ProjGnomAlgebraic:
            param = 1.0f;
            sscanf(params.c_str(), "c=%f", &param);
            _imp = new GnomonicAlgebraic<FLOAT>(param);
            _forward = GnomonicAlgebraic<FLOAT>::forward;
            _inverse = GnomonicAlgebraic<FLOAT>::inverse;
            _valid = true;
            break;
        case ProjGnomPoly5:
            param = 0.0f;
            sscanf(params.c_str(), "ii=%f", &param);
            _imp = new GnomonicPoly5<FLOAT>(param > 0.0f);
            _forward = GnomonicPoly5<FLOAT>::forward;
            _inverse = GnomonicPoly5<FLOAT>::inverse;
            _valid = true;
            break;
        case ProjGnomSmoothstep:
            param = 1.0f;
            sscanf(params.c_str(), "c=%f", &param);
            _imp = new GnomonicSmoothstep<FLOAT>(param);
            _forward = GnomonicSmoothstep<FLOAT>::forward;
            _inverse = GnomonicSmoothstep<FLOAT>::inverse;
            _valid = true;
            break;
        case ProjGnomEveritt:
            param = 1.0f;
            sscanf(params.c_str(), "c=%f", &param);
            _imp = new GnomonicEveritt<FLOAT>(param);
            _forward = GnomonicEveritt<FLOAT>::forward;
            _inverse = GnomonicEveritt<FLOAT>::inverse;
            _valid = true;
            break;
        }
    }

    bool is_valid() const { return _valid; }

    void forward(FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y) const
    {
        if (!possibly_in_cubeside(center, lat, lon)) {
            x = y = NAN;
        } else {
            _forward(this, _imp, lat, lon, x, y);
        }
    }

    void inverse(FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon) const
    {
        _inverse(this, _imp, x, y, lat, lon);
    }

    void derivatives(FLOAT lat, FLOAT lon,
            FLOAT& dx_dlat, FLOAT& dx_dlon, FLOAT& dy_dlat, FLOAT& dy_dlon, FLOAT delta = 1e-5) const
    {
        FLOAT mlat = lat - delta;
        FLOAT plat = lat + delta;
        FLOAT mlon = lon - delta;
        FLOAT plon = lon + delta;

        if (mlat < -m_pi_2<FLOAT>() || plat > m_pi_2<FLOAT>())
            goto fail;
        if (mlon < -m_pi<FLOAT>())
            mlon += 2 * m_pi<FLOAT>();
        if (plon > m_pi<FLOAT>())
            plon -= 2 * m_pi<FLOAT>();

        FLOAT x_mlat_mlon, x_mlat_plon, x_plat_mlon, x_plat_plon;
        FLOAT y_mlat_mlon, y_mlat_plon, y_plat_mlon, y_plat_plon;

        forward(mlat, mlon, x_mlat_mlon, y_mlat_mlon);
        forward(mlat, plon, x_mlat_plon, y_mlat_plon);
        forward(plat, mlon, x_plat_mlon, y_plat_mlon);
        forward(plat, plon, x_plat_plon, y_plat_plon);

        if (!isfinite(x_mlat_mlon) || !isfinite(x_mlat_plon) || !isfinite(x_plat_mlon) || !isfinite(x_plat_plon)
                || !isfinite(y_mlat_mlon) || !isfinite(y_mlat_plon) || !isfinite(y_plat_mlon) || !isfinite(y_plat_plon))
            goto fail;

        dx_dlat = (- x_plat_plon + x_mlat_plon + x_mlat_mlon - x_plat_mlon) / (4 * delta);
        dx_dlon = (+ x_plat_plon + x_mlat_plon - x_mlat_mlon - x_plat_mlon) / (4 * delta);
        dy_dlat = (+ y_plat_plon - y_mlat_plon - y_mlat_mlon + y_plat_mlon) / (4 * delta);
        dy_dlon = (- y_plat_plon - y_mlat_plon + y_mlat_mlon + y_plat_mlon) / (4 * delta);
        return;

fail:
        dx_dlat = dx_dlon = dy_dlat = dy_dlon = NAN;
    }

    void analysis(FLOAT lat, FLOAT lon,
            FLOAT& h, FLOAT& k, FLOAT& theta_prime,
            FLOAT& a, FLOAT& b,
            FLOAT& omega, FLOAT& s,
            FLOAT delta = 1e-5) const
    {
        FLOAT dx_dlat, dx_dlon, dy_dlat, dy_dlon;

        derivatives(lat, lon, dx_dlat, dx_dlon, dy_dlat, dy_dlon, delta);
        if (!isfinite(dx_dlat) || !isfinite(dx_dlon) || !isfinite(dy_dlat) || !isfinite(dy_dlon)) {
            h = k = theta_prime = a = b = omega = s = NAN;
            return;
        }

        FLOAT coslat = cos(lat);
        h = hypot(dx_dlat, dy_dlat);
        k = hypot(dx_dlon, dy_dlon) / coslat;
        FLOAT sin_theta_prime = (dy_dlat * dx_dlon - dx_dlat * dy_dlon) / (h * k * coslat);
        s = h * k * sin_theta_prime; // equivalent to s = a * b
        theta_prime = asin(min(max(sin_theta_prime, static_cast<FLOAT>(-1)), static_cast<FLOAT>(+1)));
        FLOAT a_prime = sqrt(max(h * h + k * k + 2 * s, static_cast<FLOAT>(0)));
        FLOAT b_prime = sqrt(max(h * h + k * k - 2 * s, static_cast<FLOAT>(0)));
        a = (a_prime + b_prime) / 2;
        b = (a_prime - b_prime) / 2;
        omega = 2 * asin(min(max(b_prime / a_prime, static_cast<FLOAT>(-1)), static_cast<FLOAT>(+1)));
    }
};

}

#endif
