cmake_minimum_required(VERSION 3.5)
project(_SNPLIB CXX)
set(CMAKE_CXX_STANDARD 14)
add_subdirectory(pybind11)

add_definitions(-DUSE_OPENBLAS)

add_library(_SNPLIB MODULE src/snp.cc src/data_manage.cc src/statistics.cc src/relationships.cc src/genetic_variances.cc src/gwas.cc src/uni_lmm.cc src/multi_lmm.cc)