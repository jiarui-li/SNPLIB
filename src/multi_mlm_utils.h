#ifndef SNPLIB_MULTI_MLM_UTILS_H
#define SNPLIB_MULTI_MLM_UTILS_H

#include <atomic>
#include <thread>
#include <vector>

#include "line_search.h"
#include "multi_mlm_reml.h"
#include "multi_mlm_rml.h"

extern "C" {
void CalcMMLMSigmas(const double *covariates, const double *lambda,
                    double *traits, double *res, double *vars,
                    size_t num_samples, size_t num_covariates, size_t num_dims,
                    size_t num_traits, double *sigma_e, double *sigma_g,
                    size_t num_threads);
void CalcRMMLMSigmas(const double *lambda, double *res, double *vars,
                     size_t num_samples, size_t num_covariates, size_t num_dims,
                     size_t num_traits, double *sigma_e, double *sigma_g,
                     size_t num_threads);
void CalcMMLMFull(const double *covariates, const double *lambda,
                  double *traits, double *res, double *vars, size_t num_samples,
                  size_t num_covariates, size_t num_dims, size_t num_traits,
                  double *sigma_e, double *sigma_g, double *return_res,
                  size_t num_threads);
};

#endif  // SNPLIB_MULTI_MLM_UTILS_H