cmake_minimum_required(VERSION 3.31.6)

project(myproject LANGUAGES CXX VERSION 0.0.1)
set(CMAKE_CXX_STANDARD 14)


add_executable(PhragmoplastCurvatureCalculator
        main.cpp
        main.h)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIR})
target_link_libraries(PhragmoplastCurvatureCalculator GSL::gsl GSL::gslcblas)

