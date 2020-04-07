#include "genetic_variances.h"

namespace snplib {

void CalcUniLMM(const double *traits, const double *covariates,
                const double *lambda, size_t num_samples, size_t num_covariates,
                size_t num_traits, double *vars, double *res,
                size_t num_threads) {}
void CalcMLMMSigmas(const double *covariates, const double *lambda,
                    double *traits, double *res, double *vars,
                    size_t num_samples, size_t num_covariates, size_t num_dims,
                    size_t num_traits, double *sigma_e, double *sigma_g,
                    size_t num_threads) {}
void CalcRMLMMSigmas(const double *lambda, double *res, double *vars,
                     size_t num_samples, size_t num_covariates, size_t num_dims,
                     size_t num_traits, double *sigma_e, double *sigma_g,
                     size_t num_threads) {}
}  // namespace snplib