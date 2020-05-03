#include "gwas.h"

namespace {
std::atomic_size_t index;
const std::array<double, 4> geno_table{0.0, 0.0, 1.0, 2.0};
const std::array<double, 4> mask_table{1.0, 0.0, 1.0, 1.0};
}  // namespace

namespace snplib {
void CalcLinearRegressionThread(const uint8_t *geno, size_t num_samples,
                                size_t num_snps, const double *covariates,
                                size_t num_covariates, const double *trait,
                                double *betas, double *stats) {
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + num_samples_left > 0 ? 1 : 0;
  auto *mask_d = new double[num_samples];
  auto *cov = new double[num_samples * (num_covariates + 1)];
  auto *geno_d = cov + num_samples * num_covariates;
  auto *R = new double[(num_covariates + 1) * (num_covariates + 1)];
  auto m = static_cast<int32_t>(num_samples);
  auto n = static_cast<int32_t>(num_covariates + 1);
  auto local_index = index++;
  auto *tau = new double[num_covariates];
  double w1;
  LAPACKE_dgeqrf_work(LAPACK_COL_MAJOR, m, n, nullptr, m, nullptr, &w1, -1);
  double w2;
  LAPACKE_dorgqr_work(LAPACK_COL_MAJOR, m, n, n, nullptr, m, nullptr, &w2, -1);
  auto lwork = static_cast<int32_t>(w1 > w2 ? w1 : w2);
  auto *work = new double[lwork];
  while (local_index < num_snps) {
    std::copy(covariates, covariates + num_samples * num_covariates, cov);
    UnpackGeno(geno + local_index * num_bytes, num_samples, geno_table,
               mask_table, geno_d, mask_d);
    LAPACKE_dgeqrf_work(LAPACK_COL_MAJOR, m, n, cov, m, tau, work, lwork);
    for (size_t i = 0; i < num_covariates + 1; ++i) {
      std::copy(cov + i * num_samples,
                cov + i * num_samples + num_covariates + 1,
                R + i * (num_covariates + 1));
    }
    LAPACKE_dorgqr_work(LAPACK_COL_MAJOR, m, n, n, cov, m, tau, work, lwork);
  }

  delete[] mask_d;
  delete[] cov;
  delete[] R;
  delete[] tau;
  delete[] work;
}
void CalcLinearRegressionGWAS(const uint8_t *geno, size_t num_samples,
                              size_t num_snps, const double *covariates,
                              size_t num_covariates, const double *trait,
                              double *betas, double *stats,
                              size_t num_threads) {}
void CalcLogistircGWAS(const uint8_t *geno, size_t num_samples, size_t num_snps,
                       const double *covariates, size_t num_covariates,
                       const double *trait, double *betas, double *stats,
                       size_t num_threads) {}
void CalcCCAGWAS(const uint8_t *geno, size_t num_samples, size_t num_snps,
                 const double *trait, size_t num_dims, double *betas,
                 double *stats, size_t num_threads) {}
void CalcUniLMMGWAS(const uint8_t *geno, size_t num_samples, size_t num_snps,
                    const double *covariates, size_t num_covariates,
                    const double *trait, double *betas, double *fstats,
                    double *dfs, size_t num_threads) {}
}  // namespace snplib