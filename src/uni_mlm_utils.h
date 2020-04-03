#ifndef SNPLIB_UNI_MLM_UTILS_H
#define SNPLIB_UNI_MLM_UTILS_H

#include <atomic>
#include <thread>
#include <vector>

#include "line_search.h"
#include "uni_mlm.h"

extern "C" {
void CalcUniMLM(const double *traits, const double *covariates,
                const double *lambda, size_t num_samples, size_t num_covariates,
                size_t num_traits, double *vars, double *res,
                size_t num_threads);
void CalcUniMLMGWAS(const double *trait, const double *geno_d,
                    const double *covariates, const double *lambda,
                    size_t num_samples, size_t num_covariates, size_t num_snps,
                    double *betas, double *fstats, double *dfs,
                    size_t num_threads);
};

#endif  // SNPLIB_UNI_MLM_UTILS_H