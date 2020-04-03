#ifndef SNPLIB_LOGISTIC_REGRESS_H
#define SNPLIB_LOGISTIC_REGRESS_H

#ifndef USE_OPENBLAS
#include <mkl.h>
#else
#include <cblas.h>
#include <lapacke.h>
#endif
#include <algorithm>
#include <cmath>

namespace snplib {
template <int P>
class LogisticRegress {
 private:
  const double p = static_cast<double>(P);
  size_t num_samples_;
  size_t num_covariates_;
  int32_t m;
  int32_t n;
  int32_t lwork;
  double *work;

  double *X;
  double *old_beta;
  double *beta;
  double *s;
  double *y;
  double *w;
  double *z;
  double *u;

 public:
  LogisticRegress(size_t num_samples, size_t num_covariates)
      : num_samples_(num_samples), num_covariates_(num_covariates) {
    X = new double[num_samples * num_covariates];
    old_beta = new double[num_covariates];
    beta = new double[num_covariates];
    s = new double[num_covariates];
    y = new double[num_samples];
    w = new double[num_samples];
    z = new double[num_samples];
    u = new double[num_samples];
    m = static_cast<int32_t>(num_samples);
    n = static_cast<int32_t>(num_covariates);
    double w;
    LAPACKE_dgels_work(LAPACK_COL_MAJOR, 'N', m, n, 1, nullptr, m, nullptr, m,
                       &w, -1);
    lwork = static_cast<int32_t>(w);
    work = new double[lwork];
  }
  ~LogisticRegress() noexcept {
    delete[] X;
    delete[] old_beta;
    delete[] beta;
    delete[] s;
    delete[] y;
    delete[] w;
    delete[] z;
    delete[] u;
    delete[] work;
  }
  double CalcLikelihood(const double *covariates, const double *response,
                        const double *mask) {
    cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, 1.0, covariates, m, beta, 1,
                0.0, z, 1);
    for (size_t i = 0; i < num_samples_; ++i) {
      y[i] = z[i] >= 0 ? z[i] + std::log(1.0 + std::exp(-z[i]))
                       : std::log(1.0 + std::exp(z[i]));
    }
    for (size_t i = 0; i < num_samples_; ++i) {
      u[i] = p * std::exp(z[i] - y[i]);
    }
    double result = 0.0;
    for (size_t i = 0; i < num_samples_; ++i) {
      result += (response[i] * z[i] - p * y[i]) * mask[i];
    }
    return result;
  }
  void Estimate(const double *covariates, const double *response,
                const double *mask) noexcept {
    std::fill(beta, beta + num_covariates_, 0.0);
    beta[0] = 1.0;
    double f_old = CalcLikelihood(covariates, response, mask);
    std::copy(beta, beta + num_covariates_, old_beta);
    for (size_t l = 0; l < 100 * num_covariates_; ++l) {
      for (size_t i = 0; i < num_samples_; ++i) {
        w[i] = std::sqrt(p) * std::exp(0.5 * z[i] - y[i]);
      }
      for (size_t i = 0; i < num_samples_; ++i) {
        z[i] *= w[i] * mask[i];
      }
      for (size_t i = 0; i < num_samples_; ++i) {
        z[i] += (response[i] - u[i]) * mask[i] / w[i];
      }
      for (size_t i = 0; i < num_covariates_; ++i) {
        auto *tmp_X = X + i * num_samples_;
        auto *tmp_C = covariates + i * num_samples_;
        for (size_t j = 0; j < num_samples_; ++j) {
          tmp_X[j] = tmp_C[j] * w[j] * mask[j];
        }
      }
      LAPACKE_dgels_work(LAPACK_COL_MAJOR, 'N', m, n, 1, X, m, z, m, work,
                         lwork);
      std::copy(z, z + num_covariates_, beta);
      double f_new = CalcLikelihood(covariates, response, mask);
      for (size_t i = 0; i < num_covariates_; ++i) {
        s[i] = beta[i] - old_beta[i];
      }
      while ((f_new < f_old) &&
             (std::fabs(f_new - f_old) > 1e-10 * (1 + std::fabs(f_old)))) {
        for (size_t i = 0; i < num_covariates_; ++i) {
          s[i] /= 2.0;
          beta[i] = old_beta[i] + s[i];
        }
        f_new = CalcLikelihood(covariates, response, mask);
      }
      if (std::fabs(f_new - f_old) < 1e-10 * (1 + std::fabs(f_old))) {
        for (size_t i = 0; i < num_samples_; ++i) {
          w[i] = std::sqrt(p) * std::exp(0.5 * z[i] - y[i]);
        }
        break;
      }
      f_old = f_new;
      std::copy(beta, beta + num_covariates_, old_beta);
    }
  }
  const double *GetU() const { return u; }
  const double *GetW() const { return w; }
  const double *GetBeta() const { return beta; }
};
}  // namespace snplib

#endif  // SNPLIB_LOGISTIC_REGRESS_H