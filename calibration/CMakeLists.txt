# exec src template

add_executable(preprocessing)

target_sources(preprocessing PRIVATE
    preprocessing.cc
    underlying.cc
    options.cc
    model.cc
    black_scholes.cc
    heston.cc
    utils/functions.cc
    utils/parser.cc
)

target_include_directories(preprocessing PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${EXTERNAL_DIR}/levmar
    ${HOMEBREW_DIR}/eigen/include/eigen3
    ${HOMEBREW_DIR}/boost/include
    ${HOMEBREW_DIR}/catch2/include/catch2
    ${HOMEBREW_DIR}/lapack/include
    ${HOMEBREW_DIR}/openblas/include
)

target_link_libraries(preprocessing PRIVATE
    levmar
    lapack
    openblas
    m
    boost_math_c99
)

target_compile_features(preprocessing PUBLIC cxx_std_17)

set_target_properties(preprocessing PROPERTIES CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)

target_compile_definitions(preprocessing PRIVATE EIGEN_USE_LAPACKE)