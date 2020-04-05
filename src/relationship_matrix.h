#ifndef _SNPLIB_SRC_RELATIONSHIP_MATRIX_H_
#define _SNPLIB_SRC_RELATIONSHIP_MATRIX_H_

#include <algorithm>
#include <thread>
#include <vector>
#ifdef _MSC_VER
#include <intrin.h>
#else
#include <x86intrin.h>
#endif
#include "logistic_regression.h"
#include "snp.h"

namespace snplib {
void CalcIBSMatrix(uint8_t *geno, size_t num_samples, size_t num_snps,
                   double *ibs, size_t num_threads);
void CalcGRMMatrix(uint8_t *geno, const double *af, size_t num_samples,
                   size_t num_snps, double *grm, size_t num_threads);
void CalcGCTADiagonal(uint8_t *geno, const double *af, size_t num_samples,
                      size_t num_snps, double *diagnol, size_t num_threads);
void CalcKINGMatrix(uint8_t *geno, size_t num_samples, size_t num_snps,
                    double *matrix, size_t num_threads);
void CalcUGRMMatrix(uint8_t *geno, size_t num_samples, size_t num_snps,
                    double *matrix, size_t num_threads);
void CalcAdjustedGRM(uint8_t *geno, size_t num_samples, size_t num_snps,
                     double *covariates, size_t num_covariates, double *matrix,
                     double *gcta_diag, size_t num_threads);
void CalcAdmixedGRM(uint8_t *geno, size_t num_samples, size_t num_snps,
                    double *pop_af, double *pop, size_t num_pops,
                    double *matrix, double *gcta_diag, size_t num_threads);
}  // namespace snplib

#endif  //_SNPLIB_SRC_RELATIONSHIP_MATRIX_H_
