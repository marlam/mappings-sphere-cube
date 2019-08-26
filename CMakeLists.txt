# Copyright (C) 2015, 2016, 2018, 2019
# Computer Graphics Group, University of Siegen
# Written by Martin Lambers <martin.lambers@uni-siegen.de>
#
# Copying and distribution of this file, with or without modification, are
# permitted in any medium without royalty provided the copyright notice and this
# notice are preserved. This file is offered as-is, without any warranty.

cmake_minimum_required(VERSION 2.8)
set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_STANDARD 11)
if(CMAKE_COMPILER_IS_GNUCXX)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif()

# Required libraries
find_package(Qt5Gui REQUIRED)

# The programs
include_directories(${CMAKE_SOURCE_DIR}/libspherecube)

add_executable(precision-test precision-test.cpp common.hpp)
install(TARGETS precision-test RUNTIME DESTINATION bin)

add_executable(analyzer analyzer.cpp common.hpp pngvis.hpp pngvis.cpp)
target_link_libraries(analyzer Qt5::Gui)
install(TARGETS analyzer RUNTIME DESTINATION bin)

qt5_add_resources(PLOTTER_RESOURCES plotter.qrc)
add_executable(plotter plotter.cpp common.hpp
	pdfplot.hpp pdfplot.cpp ${PLOTTER_RESOURCES})
target_link_libraries(plotter Qt5::Gui)
install(TARGETS plotter RUNTIME DESTINATION bin)
