# CMAKE minimum version required
cmake_minimum_required(VERSION 3.25)

# Intro project
project(calibration
    DESCRIPTION "C++ Heston (LSVM) Calibration"
    VERSION 0.0.1
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-Wall)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(HOMEBREW_DIR "/opt/homebrew/opt")
set(EXTERNAL_DIR "$ENV{HOME}/Documents/libs")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
find_package(Catch2 3 REQUIRED)

link_directories(
    ${EXTERNAL_DIR}/levmar/lib
    ${HOMEBREW_DIR}/boost/lib
    ${HOMEBREW_DIR}/lapack/lib
    ${HOMEBREW_DIR}/openblas/lib
    ${HOMEBREW_DIR}/catch2/lib
)

include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${EXTERNAL_DIR}/levmar
    ${HOMEBREW_DIR}/eigen/include/eigen3
    ${HOMEBREW_DIR}/boost/include
    ${HOMEBREW_DIR}/catch2/include/catch2
    ${HOMEBREW_DIR}/lapack/include
    ${HOMEBREW_DIR}/openblas/include
)

add_subdirectory(src)