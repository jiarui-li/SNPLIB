#ifndef SNPLIB_LOGISTIC_UTILS_H
#define SNPLIB_LOGISTIC_UTILS_H

#include <thread>
#include <vector>

#include "logistic_regress.h"
#include "snp.h"

extern "C" {
void UnpackResGeno(uint8_t *geno, size_t num_samples, size_t num_snps,
                   double *covariates, size_t num_covariates, double *geno_d,
                   size_t num_threads);
void CalcLogisticGWAS(const double *trait, const double *covariates,
                      uint8_t *geno, size_t num_samples, size_t num_covariates,
                      size_t num_snps, double *chi2stat, size_t num_threads);
void CalcAdjustedAF(uint8_t *geno, size_t num_samples, size_t num_snps,
                    double *covariates, size_t num_covariates, double *af,
                    size_t num_threads);
void CalcAdjustedMAF(uint8_t *geno, size_t num_samples, size_t num_snps,
                     double *covariates, size_t num_covariates, double *min_maf,
                     size_t num_threads);
void CalcSNPGFT(const double *covariates, const double *af, uint8_t *geno,
                size_t num_samples, size_t num_covariates, size_t num_snps,
                double *gft, size_t num_threads);
};
#endif  // SNPLIB_LOGISTIC_UTILS_H