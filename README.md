# mappings-sphere-cube

This is source code originally written for the paper

[M. Lambers, Survey of Cube Mapping Methods in Interactive Computer Graphics. The Visual Computer (2019). DOI 10.1007/s00371-019-01708-4](http://dx.doi.org/10.1007/s00371-019-01708-4)
([pre-print version](https://marlam.de/publications/cubemaps/))

All mapping methods are implemented in the included libspherecube library.
No external libraries are required.

The evaluation methods are implemented in the analyzer, plotter, and
precision-test utilities. The modules pngvis (for analyzer) and pdfplot (for
plotter) require Qt.

The file common.hpp defines the list of projections.
