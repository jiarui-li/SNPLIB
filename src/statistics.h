#ifndef SNPLIB_SRC_STATISTICS_H_
#define SNPLIB_SRC_STATISTICS_H_

#include <algorithm>
#include <array>
#include <atomic>
#include <cstring>
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
void CalcAlleleFrequencies(const uint8_t *geno, size_t num_samples,
                           size_t num_snps, double *af);
void CalcMissing(const uint8_t *geno, size_t num_samples, size_t num_snps,
                 double *missing);
void CalcAdjustedAF(const uint8_t *geno, size_t num_samples, size_t num_snps,
                    double *covariates, size_t num_covariates, double *af,
                    size_t num_threads);
void CalcAdjustedMAF(const uint8_t *geno, size_t num_samples, size_t num_snps,
                     double *covariates, size_t num_covariates, double *min_maf,
                     size_t num_threads);
}  // namespace snplib

#endif  // SNPLIB_SRC_STATISTICS_H_