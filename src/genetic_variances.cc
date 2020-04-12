#include "genetic_variances.h"

namespace {
std::atomic_size_t ind;
}

namespace snplib {
void CalcUniLMMThread(const double *traits, const double *covariates,
                      const double *lambda, size_t num_samples,
                      size_t num_covariates, size_t num_traits, double *vars,
                      double *res) {
  UniLMM worker(lambda, covariates, num_samples, num_covariates);
  LineSearch<UniLMM> searcher(worker);
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
void CalcUniLMM(const double *traits, const double *covariates,
                const double *lambda, size_t num_samples, size_t num_covariates,
                size_t num_traits, double *vars, double *res,
                size_t num_threads) {
  std::vector<std::thread> workers(num_threads);
  ind = 0;
  for (size_t i = 0; i < num_threads; ++i) {
    workers[i] =
        std::thread(CalcUniLMMThread, traits, covariates, lambda, num_samples,
                    num_covariates, num_traits, vars, res);
  }
  for (auto &&iter : workers) {
    iter.join();
  }
}
void CalcMLMMSigmasThread(const double *covariates, const double *lambda,
                          double *traits, double *res, double *vars,
                          size_t num_samples, size_t num_covariates,
                          size_t num_dims, size_t num_traits, double *sigma_e,
                          double *sigma_g) {
  MultiLMMREML worker(lambda, covariates, num_samples, num_covariates,
                      num_dims);
  LineSearch<MultiLMMREML> searcher(worker);
  size_t local_ind = ind++;
  while (local_ind < num_traits) {
    worker.CalcInitialGuess(traits + local_ind * num_samples * num_dims,
                            res + local_ind * num_samples * num_dims,
                            vars + local_ind * 2 * num_dims);
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
    for (size_t l = 0; l < 50 * num_dims * num_dims; ++l) {
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
    worker.GetSigmaE(sigma_e + local_ind * num_dims * num_dims);
    worker.GetSigmaG(sigma_g + local_ind * num_dims * num_dims);
    local_ind = ind++;
  }
}
void CalcMLMMSigmas(const double *covariates, const double *lambda,
                    double *traits, double *res, double *vars,
                    size_t num_samples, size_t num_covariates, size_t num_dims,
                    size_t num_traits, double *sigma_e, double *sigma_g,
                    size_t num_threads) {
  set_num_threads(1);
  std::vector<std::thread> workers(num_threads);
  ind = 0;
  for (size_t i = 0; i < num_threads; ++i) {
    workers[i] = std::thread(CalcMLMMSigmasThread, covariates, lambda, traits,
                             res, vars, num_samples, num_covariates, num_dims,
                             num_traits, sigma_e, sigma_g);
  }
  for (auto &&iter : workers) {
    iter.join();
  }
}
void CalcRMLMMSigmasThread(const double *lambda, double *res, double *vars,
                           size_t num_samples, size_t num_covariates,
                           size_t num_dims, size_t num_traits, double *sigma_e,
                           double *sigma_g) {
  MultiLMMRML worker(lambda, num_samples, num_covariates, num_dims);
  LineSearch<MultiLMMRML> searcher(worker);
  size_t local_ind = ind++;
  while (local_ind < num_traits) {
    worker.CalcInitialGuess(res + local_ind * num_dims * num_samples,
                            vars + local_ind * 2 * num_dims);
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
    for (size_t l = 0; l < 50 * num_dims * num_dims; ++l) {
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
    worker.GetSigmaE(sigma_e + local_ind * num_dims * num_dims);
    worker.GetSigmaG(sigma_g + local_ind * num_dims * num_dims);
    local_ind = ind++;
  }
}
void CalcRMLMMSigmas(const double *lambda, double *res, double *vars,
                     size_t num_samples, size_t num_covariates, size_t num_dims,
                     size_t num_traits, double *sigma_e, double *sigma_g,
                     size_t num_threads) {
  set_num_threads(1);
  std::vector<std::thread> workers(num_threads);
  ind = 0;
  for (size_t i = 0; i < num_threads; ++i) {
    workers[i] =
        std::thread(CalcRMLMMSigmasThread, lambda, res, vars, num_samples,
                    num_covariates, num_dims, num_traits, sigma_e, sigma_g);
  }
  for (auto &&iter : workers) {
    iter.join();
  }
}
}  // namespace snplib