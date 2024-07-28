# CMAKE minimum version required
cmake_minimum_required(VERSION 3.12)

# Intro project
project(calibration
    DESCRIPTION "C++ Heston (LSVM) Calibration"
    VERSION 0.0.1
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-Wall -O3 -march=native -Wno-class-memaccess)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(__INCLUDE_PATH__ "/usr/include")
set(__LIB_PATH__ "/usr/lib")
set(__LOCAL_PATH__ "/usr/local")


link_directories(
    ${__LIB_PATH__}
    ${__LIB_PATH__}/aarch64-linux-gnu

)

include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${__INCLUDE_PATH__}
    ${__INCLUDE_PATH__}/eigen3
    ${__INCLUDE_PATH__}/boost
    ${__INCLUDE_PATH__}/lapack
    ${__INCLUDE_PATH__}/openblas
)

add_subdirectory(calibration)