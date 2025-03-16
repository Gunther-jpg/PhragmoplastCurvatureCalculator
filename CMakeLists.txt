cmake_minimum_required(VERSION 3.27.8)

project(myproject LANGUAGES CXX VERSION 0.0.1)
set(CMAKE_CXX_STANDARD 23)


add_executable(PhragmoplastCurvatureCalculator
        main.cpp
        main.hpp)

find_package(GSL REQUIRED)
find_package(NLopt REQUIRED)
include_directories(${GSL_INCLUDE_DIR} ${NLopt_INCLUDE_DIR})
target_link_libraries(PhragmoplastCurvatureCalculator GSL::gsl GSL::gslcblas NLopt::nlopt)

