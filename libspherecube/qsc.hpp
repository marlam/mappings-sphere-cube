#ifndef PROJ_QSC_H
#define PROJ_QSC_H

#include "spherecube.hpp"

namespace Projection {

template <typename FLOAT>
class QSC : public Implementation<FLOAT>
{
public:
    QSC();
    static void forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y);
    static void inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* imp, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon);
};


template <typename FLOAT>
QSC<FLOAT>::QSC()
{
}

/* The diagonals on a cube face square divide the cube face into four areas.
 * The top area is the QSC area of definition, the other three areas are counted
 * counterclockwise and are handled by rotation. */
#define QSC_AREA_0 0
#define QSC_AREA_1 1
#define QSC_AREA_2 2
#define QSC_AREA_3 3

#define QSC_BASED_ON_GNOMONIC 1

template <typename FLOAT>
void QSC<FLOAT>::forward(const Context<FLOAT>* ctx, const Implementation<FLOAT>* /* imp */, FLOAT lat, FLOAT lon, FLOAT& x, FLOAT& y)
{
    FLOAT theta, phi;
    FLOAT t, mu; /* nu; */
    int area;

#if !QSC_BASED_ON_GNOMONIC
    /* This implementation is for the top cube face, other faces are handled by rotation. */
    any_face_to_top_face(ctx->center, lat, lon);

    /* Convert the input lat, lon into theta, phi as used by QSC, and
     * rotate into the area of definition. */
    phi = m_pi_2<FLOAT>() - lat;
    if (lon >= m_pi_4<FLOAT>() && lon <= m_pi_2<FLOAT>() + m_pi_4<FLOAT>()) {
        area = QSC_AREA_0;
        theta = lon - m_pi_2<FLOAT>();
    } else if (lon > m_pi_2<FLOAT>() + m_pi_4<FLOAT>() || lon <= -(m_pi_2<FLOAT>() + m_pi_4<FLOAT>())) {
        area = QSC_AREA_1;
        theta = (lon > 0 ? lon - m_pi<FLOAT>() : lon + m_pi<FLOAT>());
    } else if (lon > -(m_pi_2<FLOAT>() + m_pi_4<FLOAT>()) && lon <= -m_pi_4<FLOAT>()) {
        area = QSC_AREA_2;
        theta = lon + m_pi_2<FLOAT>();
    } else {
        area = QSC_AREA_3;
        theta = lon;
    }
#else
    /* Gnomonic */
    any_face_to_front_face(ctx->center, lat, lon);
    x = tan(lon);
    y = tan(lat) / cos(lon);
    /* Compute phi, theta as used by QSC */
    phi = acos(1 / sqrt(x * x + y * y + 1));
    theta = atan2(x, -y);
    if (theta >= m_pi_4<FLOAT>() && theta <= m_pi_2<FLOAT>() + m_pi_4<FLOAT>()) {
        area = QSC_AREA_0;
        theta = theta - m_pi_2<FLOAT>();
    } else if (theta > m_pi_2<FLOAT>() + m_pi_4<FLOAT>() || theta <= -(m_pi_2<FLOAT>() + m_pi_4<FLOAT>())) {
        area = QSC_AREA_1;
        theta = (theta > 0 ? theta - m_pi<FLOAT>() : theta + m_pi<FLOAT>());
    } else if (theta > -(m_pi_2<FLOAT>() + m_pi_4<FLOAT>()) && theta <= -m_pi_4<FLOAT>()) {
        area = QSC_AREA_2;
        theta = theta + m_pi_2<FLOAT>();
    } else {
        area = QSC_AREA_3;
        theta = theta;
    }
#endif

    /* Compute mu and nu for the area of definition.
     * For mu, see Eq. (3-21) in [OL76], but note the typos:
     * compare with Eq. (3-14). For nu, see Eq. (3-38). */
    mu = atan((12 / m_pi<FLOAT>()) * (theta + acos(sin(theta) * cos(m_pi_4<FLOAT>())) - m_pi_2<FLOAT>()));
    t = sqrt((1 - cos(phi)) / (cos(mu) * cos(mu)) / (1 - cos(atan(1 / cos(theta)))));
    /* nu = atan(t);        We don't really need nu, just t, see below. */

    /* Rotate back from the area of definition to the actual area. */
    if (area == QSC_AREA_1) {
        mu += m_pi_2<FLOAT>();
    } else if (area == QSC_AREA_2) {
        mu += m_pi<FLOAT>();
    } else if (area == QSC_AREA_3) {
        mu += m_pi_2<FLOAT>() + m_pi<FLOAT>();
    }

    /* Now compute x, y from mu and nu */
    /* t = tan(nu); */
    x = t * cos(mu);
    y = t * sin(mu);
}

template <typename FLOAT>
void QSC<FLOAT>::inverse(const Context<FLOAT>* ctx, const Implementation<FLOAT>* /* imp */, FLOAT x, FLOAT y, FLOAT& lat, FLOAT& lon)
{
    FLOAT mu, /* nu, */ cosmu, tannu;
    FLOAT tantheta, theta, cosphi, phi;
    FLOAT t;
    int area;

    /* Convert the input x, y to the mu and nu angles as used by QSC, and
     * rotate into the area of definition. */
    //nu = atan(hypot(x, y));
    tannu = hypot(x, y);
    mu = atan2(y, x);
    if (x >= 0 && x >= fabs(y)) {
        area = QSC_AREA_0;
    } else if (y >= 0 && y >= fabs(x)) {
        area = QSC_AREA_1;
        mu -= m_pi_2<FLOAT>();
    } else if (x < 0 && -x >= fabs(y)) {
        area = QSC_AREA_2;
        mu = (mu < 0 ? mu + m_pi<FLOAT>() : mu - m_pi<FLOAT>());
    } else {
        area = QSC_AREA_3;
        mu += m_pi_2<FLOAT>();
    }

    /* Compute phi and theta for the area of definition.
     * The inverse projection is not described in the original paper, but some
     * good hints can be found here (as of 2011-12-14):
     * http://fits.gsfc.nasa.gov/fitsbits/saf.93/saf.9302
     * (search for "Message-Id: <9302181759.AA25477 at fits.cv.nrao.edu>") */
    t = (m_pi<FLOAT>() / 12) * tan(mu);
    tantheta = sin(t) / (cos(t) - m_sqrt1_2<FLOAT>());
    theta = atan(tantheta);
    cosmu = cos(mu);
    //tannu = tan(nu);
    cosphi = 1 - cosmu * cosmu * tannu * tannu * (1 - cos(atan(1 / cos(theta))));
    if (cosphi < -1) {
        cosphi = -1;
    } else if (cosphi > +1) {
        cosphi = +1;
    }

#if !QSC_BASED_ON_GNOMONIC
    /* Apply the result, and rotate back from the area of definition to the actual area. */
    phi = acos(cosphi);
    lat = m_pi_2<FLOAT>() - phi;
    if (area == QSC_AREA_0) {
        lon = theta + m_pi_2<FLOAT>();
    } else if (area == QSC_AREA_1) {
        lon = (theta < 0 ? theta + m_pi<FLOAT>() : theta - m_pi<FLOAT>());
    } else if (area == QSC_AREA_2) {
        lon = theta - m_pi_2<FLOAT>();
    } else /* area == QSC_AREA_3 */ {
        lon = theta;
    }
    /* This implementation is for the top cube face, other faces are handled by rotation. */
    any_face_from_top_face(ctx->center, lat, lon);
#else
    phi = acos(cosphi);
    if (area == QSC_AREA_0) {
        theta = theta + m_pi_2<FLOAT>();
    } else if (area == QSC_AREA_1) {
        theta = (theta < 0 ? theta + m_pi<FLOAT>() : theta - m_pi<FLOAT>());
    } else if (area == QSC_AREA_2) {
        theta = theta - m_pi_2<FLOAT>();
    } else /* area == QSC_AREA_3 */ {
        theta = theta;
    }
    //y = -sin(phi) * cos(theta) / cos(phi);
    //x = sin(phi) * sin(theta) / cos(phi);
    x = tan(phi) * sin(theta);
    y = tan(phi) * (-cos(theta));
    /* Inverse Gnomonic */
    lon = atan(x);
    lat = atan(y * cos(lon));
    any_face_from_front_face(ctx->center, lat, lon);
#endif
}

}

#endif
