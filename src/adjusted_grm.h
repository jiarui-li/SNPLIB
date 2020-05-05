#ifndef SNPLIB_SRC_ADJUSTED_GRM_H_
#define SNPLIB_SRC_ADJUSTED_GRM_H_

#include <thread>
#include <vector>

#include "logistc_regression.h"
#include "snp.h"

namespace snplib {
void CalcAdjustedGRM(const uint8_t *geno, size_t num_samples, size_t num_snps,
                     double *covariates, size_t num_covariates, double *matrix,
                     double *gcta_diag, size_t num_threads);
void CalcAdmixedGRM(const uint8_t *geno, size_t num_samples, size_t num_snps,
                    double *pop_af, double *pop, size_t num_pops,
                    double *matrix, double *gcta_diag, size_t num_threads);
}  // namespace snplib

#endif  // SNPLIB_SRC_ADJUSTED_GRM_H_
