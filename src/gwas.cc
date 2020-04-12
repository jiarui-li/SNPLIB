#include "gwas.h"

namespace snplib {
std::atomic_size_t ind;
void CalcLinearGWAS(const double *trait, const double *covariates,
                    uint8_t *geno, size_t num_samples, size_t num_covariates,
                    size_t num_snps, double *chi2stat, size_t num_threads) {
  //  TODO
}
void CalcLogisticGWASThread(const double *trait, const double *covariates,
                            uint8_t *geno, size_t num_samples,
                            size_t num_covariates, size_t num_snps,
                            double *chi2stat) {
  auto *cov = new double[num_samples * (num_covariates + 1)];
  auto *geno_d = cov + num_samples * num_covariates;
  auto *mask_d = new double[num_samples];
  std::copy(covariates, covariates + num_samples * num_covariates, cov);
  SNP snp(geno, num_samples);
  LogisticRegress<1> worker_full(num_samples, num_covariates + 1);
  LogisticRegress<1> worker_reduce(num_samples, num_covariates);
  auto m = static_cast<int32_t>(num_samples);
  for (size_t l = 0; l < num_snps; ++l) {
    snp.UnpackGeno(geno_d, mask_d);
    worker_full.Estimate(cov, trait, mask_d);
    worker_reduce.Estimate(covariates, trait, mask_d);
    chi2stat[l] = worker_full.CalcLikelihood(cov, trait, mask_d) -
                  worker_reduce.CalcLikelihood(covariates, trait, mask_d);
    ++snp;
  }
  delete[] geno_d;
  delete[] mask_d;
}
void CalcLogisticGWAS(const double *trait, const double *covariates,
                      uint8_t *geno, size_t num_samples, size_t num_covariates,
                      size_t num_snps, double *chi2stat, size_t num_threads) {
  set_num_threads(1);
  std::vector<std::thread> workers(num_threads);
  auto num_snps_job = num_snps / num_threads + 1;
  auto num_snps_left = num_snps % num_threads;
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + (num_samples_left > 0 ? 1 : 0);
  for (size_t i = 0; i < num_snps_left; ++i) {
    workers[i] =
        std::thread(CalcLogisticGWASThread, trait, covariates, geno,
                    num_samples, num_covariates, num_snps_job, chi2stat);
    geno += num_bytes * num_snps_job;
    chi2stat += num_snps_job;
  }
  --num_snps_job;
  for (size_t i = num_snps_left; i < num_threads; ++i) {
    workers[i] =
        std::thread(CalcLogisticGWASThread, trait, covariates, geno,
                    num_samples, num_covariates, num_snps_job, chi2stat);
    geno += num_bytes * num_snps_job;
    chi2stat += num_snps_job;
  }
  for (auto &&iter : workers) {
    iter.join();
  }
}
void CalcUniMLMGWASThread(const double *trait, const double *geno_d,
                          const double *covariates, const double *lambda,
                          size_t num_samples, size_t num_covariates,
                          size_t num_snps, double *betas, double *fstats,
                          double *dfs) {
  auto *cov = new double[(num_covariates + 1) * num_samples];
  std::copy(covariates, covariates + num_samples * num_covariates, cov);
  auto *cov_g = cov + num_samples * num_covariates;
  UniLMM worker(lambda, cov, num_samples, num_covariates + 1);
  LineSearch<UniLMM> searcher(worker);
  size_t local_ind = ind++;
  while (local_ind < num_snps) {
    std::copy(geno_d + local_ind * num_samples,
              geno_d + (local_ind + 1) * num_samples, cov_g);
    worker.CalcInitialGuess(trait);
    double a;
    double f_old = worker.CalcLikelihood();
    double f_new = f_old;
    for (size_t i = 0; i < 15; ++i) {
      f_old = f_new;
      worker.CalcEMStep();
      f_new = worker.CalcLikelihood();
    }
    f_old = worker.CalcLikelihood();
    worker.UpdateGradients();
    f_new = f_old;
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
    f_new = worker.CalcLikelihood();
    for (size_t l = 0; l < 200; ++l) {
      worker.CalcEMStep();
      f_new = worker.CalcLikelihood();
      if ((f_old - f_new) < 1e-8 * (1.0 + std::fabs(f_old))) {
        break;
      }
      f_old = f_new;
    }
    worker.UpdateGradients();
    worker.UpdateGradients();
    auto results = worker.CalcFTest();
    betas[local_ind] = results[0];
    fstats[local_ind] = results[1];
    dfs[local_ind] = results[2];
    local_ind = ind++;
  }
  delete[] cov;
}
void CalcUniLMMGWAS(const double *trait, uint8_t *geno_d,
                    const double *covariates, const double *lambda,
                    size_t num_samples, size_t num_covariates, size_t num_snps,
                    double *betas, double *fstats, double *dfs,
                    size_t num_threads) {
  // TODO
}
}  // namespace snplib