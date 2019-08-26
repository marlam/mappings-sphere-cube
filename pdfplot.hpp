#ifndef PDFPLOT_HPP
#define PDFPLOT_HPP

#include <vector>

class QPdfWriter;
class QPainter;

/*
 * Construct a QPainter that will paint into a PDF file with the specified name
 * (the ".pdf" extension will be appended).
 */
void make_pdf_painter(const std::string& name, QPdfWriter** writer, QPainter** painter);

/*
 * Finish painting and PDF writing and destroy the associated objects.
 */
void finish_pdf_painter(QPdfWriter* writer, QPainter* painter);

/*
 * Plot a set of points, either as a multiline or as a polygon.
 */
void plot_points(QPainter* painter, bool as_line, const std::vector<std::pair<float, float>>& point_list);

/*
 * Check if a point is a valid map point; allow for some tolerance around the borders.
 */
bool point_is_valid(double x, double y, double tolerance = 0.001);

#endif
