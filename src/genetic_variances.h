#ifndef SNPLIB_SRC_GENETIC_VARIANCES_H_
#define SNPLIB_SRC_GENETIC_VARIANCES_H_

#include <atomic>
#include <thread>
#include <vector>

#include "line_search.h"
#include "multi_lmm.h"
#include "uni_lmm.h"

namespace snplib {
void CalcUniLMM(const double *traits, const double *covariates,
                const double *lambda, size_t num_samples, size_t num_covariates,
                size_t num_traits, double *vars, double *res,
                size_t num_threads);
void CalcMultiLMM_REML(const double *covariates, const double *lambda,
                       const double *traits, const double *res,
                       const double *vars, size_t num_samples,
                       size_t num_covariates, size_t num_dims,
                       size_t num_traits, double *sigma_e, double *sigma_g,
                       size_t num_threads);
void CalcMultiLMM_RML(const double *lambda, const double *res,
                      const double *vars, size_t num_samples,
                      size_t num_covariates, size_t num_dims, size_t num_traits,
                      double *sigma_e, double *sigma_g, size_t num_threads);
}  // namespace snplib

#endif  // SNPLIB_SRC_GENETIC_VARIANCES_H_
