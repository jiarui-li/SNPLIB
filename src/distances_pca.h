#ifndef _SNPLIB_SRC_DISTANCES_PCA_H_
#define _SNPLIB_SRC_DISTANCES_PCA_H_

#include <cmath>
#include <thread>
#include <vector>

#include "math_lib.h"

namespace snplib {
void CalcDistancesPCA(const double *traits, size_t num_samples,
                      size_t num_traits, size_t num_dims, size_t num_components,
                      double *scores, double *vars, double *loadings,
                      size_t num_threads);
}

#endif  //_SNPLIB_SRC_DISTANCES_PCA_H_