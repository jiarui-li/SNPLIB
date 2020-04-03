#include "multi_mlm_utils.h"

namespace {
std::atomic_size_t ind;
void CalcMMLMSigmasThread(const double *covariates, const double *lambda,
                          double *traits, double *res, double *vars,
                          size_t num_samples, size_t num_covariates,
                          size_t num_dims, size_t num_traits, double *sigma_e,
                          double *sigma_g) {
  snplib::MultiMLMREML worker(lambda, covariates, num_samples, num_covariates,
                              num_dims);
  LineSearch<snplib::MultiMLMREML> searcher(worker);
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
void CalcRMMLMSigmasThread(const double *lambda, double *res, double *vars,
                           size_t num_samples, size_t num_covariates,
                           size_t num_dims, size_t num_traits, double *sigma_e,
                           double *sigma_g) {
  snplib::MultiMLMRML worker(lambda, num_samples, num_covariates, num_dims);
  LineSearch<snplib::MultiMLMRML> searcher(worker);
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
void CalcMMLMFullThread(const double *covariates, const double *lambda,
                        double *traits, double *res, double *vars,
                        size_t num_samples, size_t num_covariates,
                        size_t num_dims, size_t num_traits, double *sigma_e,
                        double *sigma_g, double *return_res) {
  snplib::MultiMLMREML worker(lambda, covariates, num_samples, num_covariates,
                              num_dims);
  LineSearch<snplib::MultiMLMREML> searcher(worker);
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
    worker.CalcRes(return_res + local_ind * num_dims * num_samples);
    local_ind = ind++;
  }
}
}  // namespace

void CalcMMLMSigmas(const double *covariates, const double *lambda,
                    double *traits, double *res, double *vars,
                    size_t num_samples, size_t num_covariates, size_t num_dims,
                    size_t num_traits, double *sigma_e, double *sigma_g,
                    size_t num_threads) {
  std::vector<std::thread> workers(num_threads);
  ind = 0;
  for (size_t i = 0; i < num_threads; ++i) {
    workers[i] = std::thread(CalcMMLMSigmasThread, covariates, lambda, traits,
                             res, vars, num_samples, num_covariates, num_dims,
                             num_traits, sigma_e, sigma_g);
  }
  for (auto &&iter : workers) {
    iter.join();
  }
}
void CalcRMMLMSigmas(const double *lambda, double *res, double *vars,
                     size_t num_samples, size_t num_covariates, size_t num_dims,
                     size_t num_traits, double *sigma_e, double *sigma_g,
                     size_t num_threads) {
  std::vector<std::thread> workers(num_threads);
  ind = 0;
  for (size_t i = 0; i < num_threads; ++i) {
    workers[i] =
        std::thread(CalcRMMLMSigmasThread, lambda, res, vars, num_samples,
                    num_covariates, num_dims, num_traits, sigma_e, sigma_g);
  }
  for (auto &&iter : workers) {
    iter.join();
  }
}
void CalcMMLMFull(const double *covariates, const double *lambda,
                  double *traits, double *res, double *vars, size_t num_samples,
                  size_t num_covariates, size_t num_dims, size_t num_traits,
                  double *sigma_e, double *sigma_g, double *return_res,
                  size_t num_threads) {
  std::vector<std::thread> workers(num_threads);
  ind = 0;
  for (size_t i = 0; i < num_threads; ++i) {
    workers[i] = std::thread(CalcMMLMFullThread, covariates, lambda, traits,
                             res, vars, num_samples, num_covariates, num_dims,
                             num_traits, sigma_e, sigma_g, return_res);
  }
  for (auto &&iter : workers) {
    iter.join();
  }
}