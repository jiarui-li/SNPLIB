#ifndef _SNPLIB_SRC_STATISTICS_H_
#define _SNPLIB_SRC_STATISTICS_H_

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
void CalcAlleleFrequencies(uint8_t *geno, size_t num_samples, size_t num_snps,
                           double *af);
void CalcMissing(uint8_t *geno, size_t num_samples, size_t num_snps,
                 double *cr);
void CalcAdjustedAF(uint8_t *geno, size_t num_samples, size_t num_snps,
                    double *covariates, size_t num_covariates, double *af,
                    size_t num_threads);
void CalcAdjustedMAF(uint8_t *geno, size_t num_samples, size_t num_snps,
                     double *covariates, size_t num_covariates, double *min_maf,
                     size_t num_threads);
}  // namespace snplib

#endif  //_SNPLIB_SRC_STATISTICS_H_
