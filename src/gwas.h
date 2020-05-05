#ifndef SNPLIB_SRC_GWAS_H_
#define SNPLIB_SRC_GWAS_H_

#include <atomic>
#include <numeric>
#include <thread>
#include <vector>

#include "line_search.h"
#include "logistc_regression.h"
#include "uni_lmm.h"
#include "unpack_geno.h"

namespace snplib {
void CalcLinearRegressionGWAS(const uint8_t *geno, size_t num_samples,
                              size_t num_snps, const double *covariates,
                              size_t num_covariates, const double *trait,
                              double *betas, double *stats, size_t num_threads);
void CalcLogisticGWAS(const uint8_t *geno, size_t num_samples, size_t num_snps,
                      const double *covariates, size_t num_covariates,
                      const double *trait, double *betas, double *stats,
                      size_t num_threads);
void CalcCCAGWAS(const uint8_t *geno, size_t num_samples, size_t num_snps,
                 const double *trait, size_t num_dims, double *betas,
                 double *stats, size_t num_threads);
void CalcUniLMMGWAS(const uint8_t *geno, size_t num_samples, size_t num_snps,
                    const double *lambda, const double *covariates,
                    size_t num_covariates, const double *trait, double *betas,
                    double *fstats, double *dfs, size_t num_threads);
}  // namespace snplib

#endif  // SNPLIB_SRC_GWAS_H_
