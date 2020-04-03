#ifndef SNPLIB_ADJUSTED_GRM_H
#define SNPLIB_ADJUSTED_GRM_H

#include <thread>
#include <vector>

#include "logistic_regress.h"
#include "snp.h"

extern "C" {
void CalcAdjustedGRM(uint8_t *geno, size_t num_samples, size_t num_snps,
                     double *covariates, size_t num_covariates, double *matrix,
                     double *gcta_diag, size_t num_threads);
void CalcAdmixedGRM(uint8_t *geno, size_t num_samples, size_t num_snps,
                    double *pop_af, double *pop, size_t num_pops,
                    double *matrix, double *gcta_diag, size_t num_threads);
};

#endif  // SNPLIB_ADJUSTED_GRM_H
