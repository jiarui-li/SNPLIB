#ifndef SNPLIB_GRM_H
#define SNPLIB_GRM_H

#include <algorithm>
#include <cmath>
#include <thread>
#include <vector>

#include "snp.h"

extern "C" {
void CalcGRMMatrix(uint8_t *geno, const double *af, size_t num_samples,
                   size_t num_snps, double *grm, size_t num_threads);
void CalcGCTADiagonal(uint8_t *geno, const double *af, size_t num_samples,
                      size_t num_snps, double *diagnol, size_t num_threads);
void UnpackGRMGeno(uint8_t *geno, const double *af, size_t num_samples,
                   size_t num_snps, double *grm_geno);
};

#endif  // SNPLIB_GRM_H
