#ifndef PNGVIS_HPP
#define PNGVIS_HPP

#include <string>

void png_vis(const std::string& name, int w, int h, const double* data,
        bool is_based_around_1, double maxval);

#endif
