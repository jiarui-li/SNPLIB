#include "gwas.h"

namespace {
std::atomic_size_t ind;
const std::array<double, 4> geno_table{0.0, 0.0, 1.0, 2.0};
const std::array<double, 4> mask_table{1.0, 0.0, 1.0, 1.0};
}  // namespace

namespace snplib {
void CalcLinearRegressionThread(const uint8_t *geno, size_t num_samples,
                                size_t num_snps, const double *covariates,
                                size_t num_covariates, const double *trait,
                                double *betas, double *stats) {
  auto *mask_d = new double[num_samples];
  auto *y = new double[num_samples];
  auto *cov = new double[num_samples * (num_covariates + 1)];
  auto *geno_d = cov + num_samples * num_covariates;
  auto *R = new double[(num_covariates + 1) * (num_covariates + 1)];
  auto *beta = new double[num_covariates + 1];
  auto m = static_cast<int32_t>(num_samples);
  auto n = static_cast<int32_t>(num_covariates + 1);
  auto local_ind = ind++;
  auto *tau = new double[num_covariates];
  double w1;
  LAPACKE_dgeqrf_work(LAPACK_COL_MAJOR, m, n, nullptr, m, nullptr, &w1, -1);
  double w2;
  LAPACKE_dorgqr_work(LAPACK_COL_MAJOR, m, n, n, nullptr, m, nullptr, &w2, -1);
  auto lwork = static_cast<int32_t>(w1 > w2 ? w1 : w2);
  auto *work = new double[lwork];
  while (local_ind < num_snps) {
    SNP snp(geno, num_samples);
    snp += local_ind;
    snp.UnpackGeno(geno_table, mask_table, geno_d, mask_d);
    std::copy(covariates, covariates + num_samples * num_covariates, cov);
    for (size_t i = 0; i < num_samples; ++i) {
      y[i] = trait[i] * mask_d[i];
    }
    for (size_t i = 0; i < num_covariates + 1; ++i) {
      auto *tmp_c = cov + i * num_samples;
      for (size_t j = 0; j < num_samples; ++j) {
        tmp_c[j] *= mask_d[j];
      }
    }
    LAPACKE_dgeqrf_work(LAPACK_COL_MAJOR, m, n, cov, m, tau, work, lwork);
    for (size_t i = 0; i < num_covariates + 1; ++i) {
      std::copy(cov + i * num_samples,
                cov + i * num_samples + num_covariates + 1,
                R + i * (num_covariates + 1));
    }
    LAPACKE_dorgqr_work(LAPACK_COL_MAJOR, m, n, n, cov, m, tau, work, lwork);
    cblas_dgemv(CblasColMajor, CblasTrans, m, n, 1.0, cov, m, y, 1, 0.0, beta,
                1);
    double var = 0.0;
    for (size_t i = 0; i < num_samples; ++i) {
      var += y[i] * y[i];
    }
    for (size_t i = 0; i < num_covariates + 1; ++i) {
      var -= beta[i] * beta[i];
    }
    var /= (std::accumulate(mask_d, mask_d + num_samples, 0.0) -
            num_covariates - 1);
    cblas_dtrsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, n, R, n,
                beta, 1);
    betas[local_ind] = beta[num_covariates];
    stats[local_ind] =
        betas[local_ind] /
        (std::sqrt(var) / R[(num_covariates + 1) * (num_covariates + 1) - 1]);
    local_ind = ind++;
  }
  delete[] mask_d;
  delete[] y;
  delete[] cov;
  delete[] R;
  delete[] beta;
  delete[] tau;
  delete[] work;
}
void CalcLinearRegressionGWAS(const uint8_t *geno, size_t num_samples,
                              size_t num_snps, const double *covariates,
                              size_t num_covariates, const double *trait,
                              double *betas, double *stats,
                              size_t num_threads) {
  std::vector<std::thread> workers;
  set_num_threads(1);
  ind = 0;
  for (size_t i = 0; i < num_threads; ++i) {
    workers.emplace_back(CalcLinearRegressionThread, geno, num_samples,
                         num_snps, covariates, num_covariates, trait, betas,
                         stats);
  }
  for (auto &&iter : workers) {
    iter.join();
  }
}
void CalcLogisticThread(const uint8_t *geno, size_t num_samples,
                        size_t num_snps, const double *covariates,
                        size_t num_covariates, const double *trait,
                        double *betas, double *stats) {
  auto *cov = new double[num_samples * (num_covariates + 1)];
  auto *geno_d = cov + num_samples * num_covariates;
  auto *mask_d = new double[num_samples];
  std::copy(covariates, covariates + num_samples * num_covariates, cov);
  LogisticRegress<1> worker_full(num_samples, num_covariates + 1);
  LogisticRegress<1> worker_reduce(num_samples, num_covariates);
  auto local_ind = ind++;
  while (local_ind < num_snps) {
    SNP snp(geno, num_samples);
    snp += local_ind;
    snp.UnpackGeno(geno_table, mask_table, geno_d, mask_d);
    worker_full.Estimate(cov, trait, mask_d);
    worker_reduce.Estimate(covariates, trait, mask_d);
    auto *beta = worker_full.GetBeta();
    betas[local_ind] = beta[num_covariates];
    stats[local_ind] = worker_full.CalcLikelihood(cov, trait, mask_d) -
                       worker_reduce.CalcLikelihood(covariates, trait, mask_d);
    local_ind = ind++;
  }
  delete[] cov;
  delete[] mask_d;
}
void CalcLogisticGWAS(const uint8_t *geno, size_t num_samples, size_t num_snps,
                      const double *covariates, size_t num_covariates,
                      const double *trait, double *betas, double *stats,
                      size_t num_threads) {
  std::vector<std::thread> workers;
  set_num_threads(1);
  ind = 0;
  for (size_t i = 0; i < num_threads; ++i) {
    workers.emplace_back(CalcLogisticThread, geno, num_samples, num_snps,
                         covariates, num_covariates, trait, betas, stats);
  }
  for (auto &&iter : workers) {
    iter.join();
  }
}
void CalcCCAThread(const uint8_t *geno, size_t num_samples, size_t num_snps,
                   const double *U, const double *S, const double *VT,
                   size_t num_dims, double *betas, double *rho2) {
  auto *geno_d = new double[num_samples];
  auto *tmp = new double[num_dims];
  auto m = static_cast<int32_t>(num_samples);
  auto n = static_cast<int32_t>(num_dims);
  auto local_ind = ind++;
  while (local_ind < num_snps) {
    SNP snp(geno, num_samples);
    snp += local_ind;
    snp.UnpackGeno(geno_table, geno_d);
    auto mu = std::accumulate(geno_d, geno_d + num_samples, 0.0);
    mu /= num_samples;
    for (size_t i = 0; i < num_samples; ++i) {
      geno_d[i] -= mu;
    }
    cblas_dgemv(CblasColMajor, CblasTrans, m, n, 1.0, U, m, geno_d, 1, 0.0, tmp,
                1);
    double a = 0.0;
    for (size_t i = 0; i < num_dims; ++i) {
      a += tmp[i] * tmp[i];
    }
    double b = 0.0;
    for (size_t i = 0; i < num_samples; ++i) {
      b += geno_d[i] * geno_d[i];
    }
    rho2[local_ind] = a / b;
    for (size_t i = 0; i < num_dims; ++i) {
      tmp[i] /= S[i];
    }
    auto *tmp_beta = betas + local_ind * num_dims;
    cblas_dgemv(CblasColMajor, CblasTrans, n, n, 1.0, VT, n, tmp, 1, 0.0,
                tmp_beta, 1);
    local_ind = ind++;
  }
  delete[] geno_d;
}
void CalcCCAGWAS(const uint8_t *geno, size_t num_samples, size_t num_snps,
                 const double *trait, size_t num_dims, double *betas,
                 double *rho2, size_t num_threads) {
  std::vector<std::thread> workers;
  set_num_threads(num_threads);
  double *U = new double[num_samples * num_dims];
  double *S = new double[num_dims];
  double *VT = new double[num_dims * num_dims];
  std::copy(trait, trait + num_samples * num_dims, U);
  auto m = static_cast<int32_t>(num_samples);
  auto n = static_cast<int32_t>(num_dims);
  LAPACKE_dgesdd(LAPACK_COL_MAJOR, 'O', m, n, U, m, S, nullptr, m, VT, n);
  set_num_threads(1);
  ind = 0;
  for (size_t i = 0; i < num_threads; ++i) {
    workers.emplace_back(CalcCCAThread, geno, num_samples, num_snps, U, S, VT,
                         num_dims, betas, rho2);
  }
  for (auto &&iter : workers) {
    iter.join();
  }
  delete[] U;
  delete[] S;
  delete[] VT;
}
void CalcCCAXThread(const uint8_t *geno, size_t num_samples, size_t num_snps,
                    const double *U, const double *S, const double *VT,
                    const double *sex, size_t num_dims, double *betas,
                    double *rho2) {
  auto *geno_d = new double[num_samples];
  auto *tmp = new double[num_dims];
  auto m = static_cast<int32_t>(num_samples);
  auto n = static_cast<int32_t>(num_dims);
  auto local_ind = ind++;
  while (local_ind < num_snps) {
    SNP snp(geno, num_samples);
    snp += local_ind;
    snp.UnpackGeno(geno_table, geno_d);
    for (size_t i = 0; i < num_samples; ++i) {
      geno_d[i] /= sex[i];
    }
    auto mu = std::accumulate(geno_d, geno_d + num_samples, 0.0);
    mu /= num_samples;
    for (size_t i = 0; i < num_samples; ++i) {
      geno_d[i] -= mu;
    }
    cblas_dgemv(CblasColMajor, CblasTrans, m, n, 1.0, U, m, geno_d, 1, 0.0, tmp,
                1);
    double a = 0.0;
    for (size_t i = 0; i < num_dims; ++i) {
      a += tmp[i] * tmp[i];
    }
    double b = 0.0;
    for (size_t i = 0; i < num_samples; ++i) {
      b += geno_d[i] * geno_d[i];
    }
    rho2[local_ind] = a / b;
    for (size_t i = 0; i < num_dims; ++i) {
      tmp[i] /= S[i];
    }
    auto *tmp_beta = betas + local_ind * num_dims;
    cblas_dgemv(CblasColMajor, CblasTrans, n, n, 1.0, VT, n, tmp, 1, 0.0,
                tmp_beta, 1);
    local_ind = ind++;
  }
  delete[] geno_d;
}
void CalcCCAGWASX(const uint8_t *geno, size_t num_samples, size_t num_snps,
                  const double *trait, const double *sex, size_t num_dims,
                  double *betas, double *rho2, size_t num_threads) {
  std::vector<std::thread> workers;
  set_num_threads(num_threads);
  double *U = new double[num_samples * num_dims];
  double *S = new double[num_dims];
  double *VT = new double[num_dims * num_dims];
  std::copy(trait, trait + num_samples * num_dims, U);
  auto m = static_cast<int32_t>(num_samples);
  auto n = static_cast<int32_t>(num_dims);
  LAPACKE_dgesdd(LAPACK_COL_MAJOR, 'O', m, n, U, m, S, nullptr, m, VT, n);
  set_num_threads(1);
  ind = 0;
  for (size_t i = 0; i < num_threads; ++i) {
    workers.emplace_back(CalcCCAXThread, geno, num_samples, num_snps, U, S, VT,
                         sex, num_dims, betas, rho2);
  }
  for (auto &&iter : workers) {
    iter.join();
  }
  delete[] U;
  delete[] S;
  delete[] VT;
}
void CCAReplicationThread(const uint8_t *geno, size_t num_samples,
                          size_t num_snps, const double *scores,
                          const double *betas, size_t num_dims, double *rho) {
  auto *geno_d = new double[num_samples];
  auto *trait = new double[num_samples];
  auto *tmp = new double[num_dims];
  auto m = static_cast<int32_t>(num_samples);
  auto n = static_cast<int32_t>(num_dims);
  auto local_ind = ind++;
  while (local_ind < num_snps) {
    SNP snp(geno, num_samples);
    snp += local_ind;
    snp.UnpackGeno(geno_table, geno_d);
    auto mu = std::accumulate(geno_d, geno_d + num_samples, 0.0);
    mu /= num_samples;
    for (size_t i = 0; i < num_samples; ++i) {
      geno_d[i] -= mu;
    }
    auto *beta = betas + local_ind * num_dims;
    cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, 1.0, scores, m, beta, 1, 0.0,
                trait, 1);
    mu = std::accumulate(trait, trait + num_samples, 0.0);
    mu /= num_samples;
    double cov = 0.0;
    double var_g = 0.0;
    double var_t = 0.0;
    for (size_t i = 0; i < num_samples; ++i) {
      trait[i] -= mu;
      cov += geno_d[i] * trait[i];
      var_g += geno_d[i] * geno_d[i];
      var_t += trait[i] * trait[i];
    }
    rho[local_ind] = cov / std::sqrt(var_g * var_t);
    local_ind = ind++;
  }
  delete[] geno_d;
  delete[] trait;
}
void CalcCCAReplication(const uint8_t *geno, size_t num_samples,
                        size_t num_snps, const double *scores,
                        const double *betas, size_t num_dims, double *rho,
                        size_t num_threads) {
  std::vector<std::thread> workers;
  set_num_threads(1);
  ind = 0;
  for (size_t i = 0; i < num_threads; ++i) {
    workers.emplace_back(CCAReplicationThread, geno, num_samples, num_snps,
                         scores, betas, num_dims, rho);
  }
  for (auto &&iter : workers) {
    iter.join();
  }
}
void CCAReplicationXThread(const uint8_t *geno, size_t num_samples,
                           size_t num_snps, const double *scores,
                           const double *betas, const double *sex,
                           size_t num_dims, double *rho) {
  auto *geno_d = new double[num_samples];
  auto *trait = new double[num_samples];
  auto *tmp = new double[num_dims];
  auto m = static_cast<int32_t>(num_samples);
  auto n = static_cast<int32_t>(num_dims);
  auto local_ind = ind++;
  while (local_ind < num_snps) {
    SNP snp(geno, num_samples);
    snp += local_ind;
    snp.UnpackGeno(geno_table, geno_d);
    for (size_t i = 0; i < num_samples; ++i) {
      geno_d[i] /= sex[i];
    }
    auto mu = std::accumulate(geno_d, geno_d + num_samples, 0.0);
    mu /= num_samples;
    for (size_t i = 0; i < num_samples; ++i) {
      geno_d[i] -= mu;
    }
    auto *beta = betas + local_ind * num_dims;
    cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, 1.0, scores, m, beta, 1, 0.0,
                trait, 1);
    mu = std::accumulate(trait, trait + num_samples, 0.0);
    mu /= num_samples;
    double cov = 0.0;
    double var_g = 0.0;
    double var_t = 0.0;
    for (size_t i = 0; i < num_samples; ++i) {
      trait[i] -= mu;
      cov += geno_d[i] * trait[i];
      var_g += geno_d[i] * geno_d[i];
      var_t += trait[i] * trait[i];
    }
    rho[local_ind] = cov / std::sqrt(var_g * var_t);
    local_ind = ind++;
  }
  delete[] geno_d;
  delete[] trait;
}
void CalcCCAReplicationX(const uint8_t *geno, size_t num_samples,
                         size_t num_snps, const double *scores,
                         const double *betas, const double *sex,
                         size_t num_dims, double *rho, size_t num_threads) {
  std::vector<std::thread> workers;
  set_num_threads(1);
  ind = 0;
  for (size_t i = 0; i < num_threads; ++i) {
    workers.emplace_back(CCAReplicationXThread, geno, num_samples, num_snps,
                         scores, betas, sex, num_dims, rho);
  }
  for (auto &&iter : workers) {
    iter.join();
  }
}
void CalcUniLMMThread(const uint8_t *geno, size_t num_samples, size_t num_snps,
                      const double *lambda, const double *V,
                      const double *covariates, size_t num_covariates,
                      const double *trait, double *betas, double *fstats,
                      double *dfs) {
  auto *cov = new double[(num_covariates + 1) * num_samples];
  std::copy(covariates, covariates + num_samples * num_covariates, cov);
  UniLMM worker(lambda, cov, num_samples, num_covariates + 1);
  LineSearch<UniLMM> searcher(worker);
  auto m = static_cast<int32_t>(num_samples);
  auto *geno_d = cov + num_samples * num_covariates;
  auto *mask_d = new double[num_samples];
  auto *y = new double[num_samples];
  auto local_ind = ind++;
  while (local_ind < num_snps) {
    SNP snp(geno, num_samples);
    snp += local_ind;
    snp.UnpackGeno(geno_table, mask_table, geno_d, mask_d);
    cblas_dgemv(CblasColMajor, CblasTrans, m, m, 1.0, V, m, y, 1, 0.0, geno_d,
                1);
    for (size_t i = 0; i < num_samples; ++i) {
      y[i] = trait[i] * mask_d[i];
    }
    for (size_t i = 0; i < num_covariates + 1; ++i) {
      auto *tmp_c = cov + i * num_samples;
      for (size_t j = 0; j < num_samples; ++j) {
        tmp_c[j] *= mask_d[j];
      }
    }
    worker.CalcInitialGuess(y);
    double a;
    double f_old = worker.CalcLikelihood();
    for (size_t i = 0; i < 15; ++i) {
      worker.CalcEMStep();
      f_old = worker.CalcLikelihood();
    }
    worker.UpdateGradients();
    double f_new = f_old;
    for (size_t l = 0; l < 200; ++l) {
      f_old = f_new;
      worker.BackupVars();
      worker.UpdateHessian();
      worker.CalcAIStep();
      a = searcher.Search(f_old);
      f_new = searcher.GetFNew();
      if ((f_old - f_new) < 1e-8 * (1.0 + std::fabs(f_old))) {
        break;
      }
    }
    auto result = worker.CalcFTest();
    betas[local_ind] = result[0];
    fstats[local_ind] = result[1];
    dfs[local_ind] = result[2];
    local_ind = ind++;
  }
  delete[] cov;
  delete[] mask_d;
  delete[] y;
}
void CalcUniLMMGWAS(const uint8_t *geno, size_t num_samples, size_t num_snps,
                    const double *lambda, const double *V,
                    const double *covariates, size_t num_covariates,
                    const double *trait, double *betas, double *fstats,
                    double *dfs, size_t num_threads) {
  std::vector<std::thread> workers;
  set_num_threads(1);
  ind = 0;
  for (size_t i = 0; i < num_threads; ++i) {
    workers.emplace_back(CalcUniLMMThread, geno, num_samples, num_snps, lambda,
                         V, covariates, num_covariates, trait, betas, fstats,
                         dfs);
  }
  for (auto &&iter : workers) {
    iter.join();
  }
}
}  // namespace snplib