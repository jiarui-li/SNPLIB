#ifndef _SNPLIB_SRC_ANCESTRY_ESTIMATION_H_
#define _SNPLIB_SRC_ANCESTRY_ESTIMATION_H_

#include <algorithm>
#include <cmath>
#include <random>

#include "math_lib.h"
#include "relationships.h"
#include "statistics.h"

namespace snplib {
void CalcPCALoadingsExact(uint8_t *geno, size_t num_samples, size_t num_snps,
                          size_t num_components, double *loadings,
                          size_t num_threads);
void CalcPCALoadingsApprox(uint8_t *geno, size_t num_samples, size_t num_snps,
                           size_t num_components, double *loadings,
                           size_t num_threads);
void ProjectPCA(uint8_t *dest_geno, size_t num_samples, size_t num_snps,
                double *loadings, size_t num_components, double *scores,
                size_t num_threads);
void CalcSUGIBSLoadingsExact(uint8_t *geno, size_t num_samples, size_t num_snps,
                             size_t num_components, double *loadings,
                             size_t num_threads);
void CalcSUGIBSLoadingsApprox(uint8_t *geno, size_t num_samples,
                              size_t num_snps, size_t num_components,
                              double *loadings, size_t num_threads);
void ProjectSUGIBS(uint8_t *src_geno, size_t num_src_samples,
                   uint8_t *dest_geno, size_t num_dest_samples, size_t num_snps,
                   double *loadings, size_t num_components, double *scores,
                   size_t num_threads);
void CalcUPCALoadingsExact(uint8_t *geno, size_t num_samples, size_t num_snps,
                           size_t num_components, double *loadings,
                           size_t num_threads);
void CalcUPCALoadingsApprox(uint8_t *geno, size_t num_samples, size_t num_snps,
                            size_t num_components, double *loadings,
                            size_t num_threads);
void ProjectUPCA(uint8_t *dest_geno, size_t num_samples, size_t num_snps,
                 double *loadings, size_t num_components, double *scores,
                 size_t num_threads);
}  // namespace snplib

#endif  // _SNPLIB_SRC_ANCESTRY_ESTIMATION_H_