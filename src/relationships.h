#ifndef _SNPLIB_SRC_RELATIONSHIPS_H_
#define _SNPLIB_SRC_RELATIONSHIPS_H_

#include <thread>
#include <vector>
#ifdef _MSC_VER
#include <intrin.h>
#else
#include <x86intrin.h>
#endif

#include "logistc_regression.h"
#include "snp.h"

namespace snplib {
void CalcAdjustedGRM(uint8_t *geno, size_t num_samples, size_t num_snps,
                     double *covariates, size_t num_covariates, double *matrix,
                     double *gcta_diag, size_t num_threads);
void CalcAdmixedGRM(uint8_t *geno, size_t num_samples, size_t num_snps,
                    double *pop_af, double *pop, size_t num_pops,
                    double *matrix, double *gcta_diag, size_t num_threads);
void CalcGRMMatrix(uint8_t *geno, const double *af, size_t num_samples,
                   size_t num_snps, double *grm, size_t num_threads);
void CalcGCTADiagonal(uint8_t *geno, const double *af, size_t num_samples,
                      size_t num_snps, double *diagnol, size_t num_threads);
void CalcIBSMatrix(uint8_t *geno, size_t num_samples, size_t num_snps,
                   double *ibs, size_t num_threads);
void CalcIBSConnection(uint8_t *src_geno, size_t num_src_samples,
                       uint8_t *dest_geno, size_t num_dest_samples,
                       size_t num_snps, double *connection, size_t num_threads);
void CalcKINGMatrix(uint8_t *geno, size_t num_samples, size_t num_snps,
                    double *matrix, size_t num_threads);
void CalcUGRMMatrix(uint8_t *geno, size_t num_samples, size_t num_snps,
                    double *matrix, size_t num_threads);
}  // namespace snplib

#endif  //_SNPLIB_SRC_RELATIONSHIPS_H_
