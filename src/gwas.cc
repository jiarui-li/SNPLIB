#include "gwas.h"

namespace snplib {

void CalcLinearGWAS(const double *trait, const double *covariates,
                    uint8_t *geno, size_t num_samples, size_t num_covariates,
                    size_t num_snps, double *chi2stat, size_t num_threads) {}
void CalcLogisticGWAS(const double *trait, const double *covariates,
                      uint8_t *geno, size_t num_samples, size_t num_covariates,
                      size_t num_snps, double *chi2stat, size_t num_threads) {}
void CalcUniLMMGWAS(const double *trait, uint8_t *geno_d,
                    const double *covariates, const double *lambda,
                    size_t num_samples, size_t num_covariates, size_t num_snps,
                    double *betas, double *fstats, double *dfs,
                    size_t num_threads) {}
}  // namespace snplib