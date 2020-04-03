#include "uni_mlm_utils.h"

namespace {
std::atomic_size_t ind;
void CalcUniMLMThread(const double *traits, const double *covariates,
                      const double *lambda, size_t num_samples,
                      size_t num_covariates, size_t num_traits, double *vars,
                      double *res) {
  snplib::UniMLM worker(lambda, covariates, num_samples, num_covariates);
  LineSearch<snplib::UniMLM> searcher(worker);
  size_t local_ind = ind++;
  while (local_ind < num_traits) {
    worker.CalcInitialGuess(traits + local_ind * num_samples);
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
    auto var = worker.GetVars();
    vars[2 * local_ind] = var[0];
    vars[2 * local_ind + 1] = var[1];
    worker.CalcRes(res + local_ind * num_samples);
    local_ind = ind++;
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
  snplib::UniMLM worker(lambda, cov, num_samples, num_covariates + 1);
  LineSearch<snplib::UniMLM> searcher(worker);
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
}  // namespace

void CalcUniMLM(const double *traits, const double *covariates,
                const double *lambda, size_t num_samples, size_t num_covariates,
                size_t num_traits, double *vars, double *res,
                size_t num_threads) {
  std::vector<std::thread> workers(num_threads);
  ind = 0;
  for (size_t i = 0; i < num_threads; ++i) {
    workers[i] =
        std::thread(CalcUniMLMThread, traits, covariates, lambda, num_samples,
                    num_covariates, num_traits, vars, res);
  }
  for (auto &&iter : workers) {
    iter.join();
  }
}
void CalcUniMLMGWAS(const double *trait, const double *geno_d,
                    const double *covariates, const double *lambda,
                    size_t num_samples, size_t num_covariates, size_t num_snps,
                    double *betas, double *fstats, double *dfs,
                    size_t num_threads) {
  std::vector<std::thread> workers(num_threads);
  ind = 0;
  for (size_t i = 0; i < num_threads; ++i) {
    workers[i] =
        std::thread(CalcUniMLMGWASThread, trait, geno_d, covariates, lambda,
                    num_samples, num_covariates, num_snps, betas, fstats, dfs);
  }
  for (auto &&iter : workers) {
    iter.join();
  }
}