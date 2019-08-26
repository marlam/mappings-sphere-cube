#include <cmath>

#include <QPdfWriter>
#include <QPainter>

#include "pdfplot.hpp"

void make_pdf_painter(const std::string& name, QPdfWriter** writer, QPainter** painter)
{
    *writer = new QPdfWriter((name + ".pdf").c_str());
    (*writer)->setPageLayout(QPageLayout(
                QPageSize(QSizeF(210.0f, 210.0f), QPageSize::Millimeter, "A4 square", QPageSize::ExactMatch),
                QPageLayout::Portrait, QMarginsF(0.0f, 0.0f, 0.0f, 0.0f), QPageLayout::Millimeter));
    (*writer)->setTitle(name.c_str());
    *painter = new QPainter(*writer);
    (*painter)->setClipRegion(QRegion(0, 0, (*painter)->device()->width(), (*painter)->device()->height(), QRegion::Rectangle));
}

void finish_pdf_painter(QPdfWriter* writer, QPainter* painter)
{
    painter->end();
    delete painter;
    delete writer;
}

void plot_points(QPainter* painter, bool as_line, const std::vector<std::pair<float, float>>& point_list)
{
    float w = painter->device()->width();
    float h = painter->device()->height();
    std::vector<QPointF> qpoints;
    for (size_t i = 0; i < point_list.size(); i++) {
        float x = point_list[i].first;
        float y = point_list[i].second;
        qpoints.push_back(QPointF((x + 1.0f) / 2.0f * w, (-y + 1.0f) / 2.0f * h));
    }
    if (as_line)
        painter->drawPolyline(qpoints.data(), qpoints.size());
    else
        painter->drawPolygon(qpoints.data(), qpoints.size());
}

bool point_is_valid(double x, double y, double tolerance)
{
    return (std::isfinite(x) && std::isfinite(y) && std::fabs(x) < 1.0 + tolerance && std::fabs(y) < 1.0 + tolerance);
}
