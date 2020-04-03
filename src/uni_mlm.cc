#include "uni_mlm.h"

namespace snplib {
double UniMLM::CalcLikelihood() {
  double v_e = vars_[0] * vars_[0];
  double v_g = vars_[1] * vars_[1];
  for (size_t i = 0; i < num_samples_; ++i) {
    h[i] = v_e + lambda_[i] * v_g;
    t[i] = trait_[i] / h[i];
    p[i] = 1.0 / h[i];
  }
  double f = 0.0;
  for (size_t i = 0; i < num_samples_; ++i) {
    f += std::log(h[i]);
    h[i] = std::sqrt(h[i]);
  }
  for (size_t i = 0; i < num_covariates_; ++i) {
    auto *tmp_q = Q + i * num_samples_;
    auto *tmp_c = covariates_ + i * num_samples_;
    for (size_t j = 0; j < num_samples_; ++j) {
      tmp_q[j] = tmp_c[j] / h[j];
    }
  }
  LAPACKE_dgeqrf_work(LAPACK_COL_MAJOR, m, n, Q, m, tau, work, lwork);
  for (size_t i = 0; i < num_covariates_; ++i) {
    f += std::log(Q[i * num_samples_ + i] * Q[i * num_samples_ + i]);
  }
  LAPACKE_dorgqr_work(LAPACK_COL_MAJOR, m, n, n, Q, m, tau, work, lwork);
  for (size_t i = 0; i < num_covariates_; ++i) {
    auto *tmp_q = Q + i * num_samples_;
    for (size_t j = 0; j < num_samples_; ++j) {
      tmp_q[j] /= h[j];
    }
  }
  cblas_dgemv(CblasColMajor, CblasTrans, m, n, 1.0, Q, m, trait_, 1, 0.0, r, 1);
  cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, -1.0, Q, m, r, 1, 1.0, t, 1);
  for (int i = 0; i < num_samples_; ++i) {
    f += trait_[i] * t[i];
  }
  return f / 2.0;
}
void UniMLM::UpdateGradients() {
  for (size_t i = 0; i < num_samples_; ++i) {
    double pp = 0.0;
    for (size_t j = 0; j < num_covariates_; ++j) {
      pp += Q[j * num_samples_ + i] * Q[j * num_samples_ + i];
    }
    p[i] -= pp;
  }
  double g_e = 0.0;
  double g_g = 0.0;
  for (size_t i = 0; i < num_samples_; ++i) {
    g_e += p[i] - t[i] * t[i];
    g_g += lambda_[i] * (p[i] - t[i] * t[i]);
  }
  gradients_[0] = vars_[0] * g_e;
  gradients_[1] = vars_[1] * g_g;
}
void UniMLM::UpdateHessian() {
  double h_ee = 0.0;
  double h_ge = 0.0;
  double h_gg = 0.0;
  for (size_t i = 0; i < num_samples_; ++i) {
    h_ee += t[i] * t[i] / h[i] / h[i];
    h_ge += lambda_[i] * t[i] * t[i] / h[i] / h[i];
    h_gg += lambda_[i] * lambda_[i] * t[i] * t[i] / h[i] / h[i];
  }
  cblas_dgemv(CblasColMajor, CblasTrans, m, n, 1.0, Q, m, t, 1, 0.0, r, 1);
  for (size_t i = 0; i < num_samples_; ++i) {
    t[i] *= lambda_[i];
  }
  cblas_dgemv(CblasColMajor, CblasTrans, m, n, 1.0, Q, m, t, 1, 0.0, h, 1);
  for (size_t i = 0; i < num_covariates_; ++i) {
    h_ee += r[i] * r[i];
    h_ge += h[i] * r[i];
    h_gg += h[i] * h[i];
  }
  hessian_[0] = 2.0 * vars_[0] * vars_[0] * h_ee;
  hessian_[1] = 2.0 * vars_[1] * vars_[0] * h_ge;
  hessian_[2] = 2.0 * vars_[1] * vars_[1] * h_gg;
}

void UniMLM::CalcAIStep() {
  double hdeter = hessian_[0] * hessian_[2] - hessian_[1] * hessian_[1];
  step_[0] = hessian_[1] * gradients_[1] - hessian_[2] * gradients_[0];
  step_[1] = hessian_[1] * gradients_[0] - hessian_[0] * gradients_[1];
  step_[0] /= hdeter;
  step_[1] /= hdeter;
}

void UniMLM::CalcEMStep() {
  double v_e = vars_[0] * vars_[0];
  double v_g = vars_[1] * vars_[1];
  for (size_t i = 0; i < num_samples_; ++i) {
    double pp = 0.0;
    for (size_t j = 0; j < num_covariates_; ++j) {
      pp += Q[j * num_samples_ + i] * Q[j * num_samples_ + i];
    }
    p[i] -= pp;
  }
  double g_e = 0.0;
  double g_g = 0.0;
  for (size_t i = 0; i < num_samples_; ++i) {
    g_e += p[i] - t[i] * t[i];
    g_g += lambda_[i] * (p[i] - t[i] * t[i]);
  }
  g_e *= v_e * v_e / m;
  g_g *= v_g * v_g / m;
  v_e -= g_e;
  v_g -= g_g;
  vars_[0] = std::sqrt(v_e);
  vars_[1] = std::sqrt(v_g);
}

UniMLM::UniMLM(const double *lambda, const double *covariates,
               size_t num_samples, size_t num_covariates)
    : num_samples_(num_samples),
      num_covariates_(num_covariates),
      lambda_(lambda),
      covariates_(covariates) {
  m = static_cast<int32_t>(num_samples);
  n = static_cast<int32_t>(num_covariates);
  Q = new double[num_covariates * num_samples];
  h = new double[num_samples];
  p = new double[num_samples];
  t = new double[num_samples];
  r = new double[num_covariates];
  tau = new double[num_covariates];
  double w1;
  LAPACKE_dgeqrf_work(LAPACK_COL_MAJOR, m, n, nullptr, m, nullptr, &w1, -1);
  double w2;
  LAPACKE_dorgqr_work(LAPACK_COL_MAJOR, m, n, n, nullptr, m, nullptr, &w2, -1);
  lwork = static_cast<int32_t>(w1 > w2 ? w1 : w2);
  work = new double[lwork];
}
UniMLM::~UniMLM() noexcept {
  delete[] Q;
  delete[] h;
  delete[] p;
  delete[] t;
  delete[] r;
  delete[] tau;
  delete[] work;
}
double UniMLM::CalcLineGradient() {
  return gradients_[0] * step_[0] + gradients_[1] * step_[1];
}
void UniMLM::CalcInitialGuess(const double *trait) {
  trait_ = trait;
  std::copy(covariates_, covariates_ + num_samples_ * num_covariates_, Q);
  LAPACKE_dgeqrf_work(LAPACK_COL_MAJOR, m, n, Q, m, tau, work, lwork);
  LAPACKE_dorgqr_work(LAPACK_COL_MAJOR, m, n, n, Q, m, tau, work, lwork);
  cblas_dgemv(CblasColMajor, CblasTrans, m, n, 1.0, Q, m, trait, 1, 0.0, r, 1);
  double var = 0.0;
  for (size_t i = 0; i < num_samples_; ++i) {
    var += trait[i] * trait[i];
  }
  for (size_t i = 0; i < num_covariates_; ++i) {
    var -= r[i] * r[i];
  }
  var /= (num_samples_ - num_covariates_);
  vars_[0] = std::sqrt(var / 2.0);
  vars_[1] = std::sqrt(var / 2.0);
}
void UniMLM::BackupVars() {
  old_vars_[0] = vars_[0];
  old_vars_[1] = vars_[1];
}
void UniMLM::UpdateVars(double a) {
  vars_[0] = old_vars_[0] + a * step_[0];
  vars_[1] = old_vars_[1] + a * step_[1];
}
void UniMLM::CalcRes(double *res) {
  for (size_t i = 0; i < num_samples_; ++i) {
    res[i] = t[i] * h[i] * h[i];
  }
}
void UniMLM::CalcBeta(double *beta) {
  double v_e = vars_[0] * vars_[0];
  double v_g = vars_[1] * vars_[1];
  for (size_t i = 0; i < num_samples_; ++i) {
    h[i] = v_e + lambda_[i] * v_g;
    t[i] = trait_[i] / h[i];
  }
  for (size_t i = 0; i < num_covariates_; ++i) {
    auto *tmp_q = Q + i * num_samples_;
    auto *tmp_c = covariates_ + i * num_samples_;
    for (size_t j = 0; j < num_samples_; ++j) {
      tmp_q[j] = tmp_c[j] / h[j];
    }
  }
  LAPACKE_dgeqrf_work(LAPACK_COL_MAJOR, m, n, Q, m, tau, work, lwork);
  LAPACKE_dormqr_work(LAPACK_COL_MAJOR, 'L', 'T', m, 1, n, Q, m, tau, t, m,
                      work, lwork);
  LAPACKE_dtrtrs(LAPACK_COL_MAJOR, 'U', 'N', 'N', n, 1, Q, m, t, m);
  std::copy(t, t + num_covariates_, beta);
}
std::array<double, 3> UniMLM::CalcFTest() {
  double v_e = vars_[0] * vars_[0];
  double v_g = vars_[1] * vars_[1];
  for (size_t i = 0; i < num_samples_; ++i) {
    h[i] = v_e + lambda_[i] * v_g;
    h[i] = std::sqrt(h[i]);
    t[i] = trait_[i] / h[i];
  }
  for (size_t i = 0; i < num_covariates_; ++i) {
    auto *tmp_q = Q + i * num_samples_;
    auto *tmp_c = covariates_ + i * num_samples_;
    for (size_t j = 0; j < num_samples_; ++j) {
      tmp_q[j] = tmp_c[j] / h[j];
    }
  }
  LAPACKE_dgeqrf_work(LAPACK_COL_MAJOR, m, n, Q, m, tau, work, lwork);
  double rr = Q[(num_covariates_ - 1) * num_samples_ + num_covariates_ - 1];
  LAPACKE_dormqr_work(LAPACK_COL_MAJOR, 'L', 'T', m, 1, n, Q, m, tau, t, m,
                      work, lwork);
  LAPACKE_dtrtrs(LAPACK_COL_MAJOR, 'U', 'N', 'N', n, 1, Q, m, t, m);
  double beta = t[num_covariates_ - 1];
  double fstat = beta * beta * rr * rr;
  LAPACKE_dorgqr_work(LAPACK_COL_MAJOR, m, n, n, Q, m, tau, work, lwork);
  auto *q_n = Q + (num_covariates_ - 1) * num_samples_;
  double g_e = 0.0;
  for (size_t i = 0; i < num_samples_; ++i) {
    g_e += q_n[i] * q_n[i] / h[i] / h[i];
  }
  g_e /= rr * rr;
  double g_g = 0.0;
  for (size_t i = 0; i < num_samples_; ++i) {
    g_g += lambda_[i] * q_n[i] * q_n[i] / h[i] / h[i];
  }
  g_g /= rr * rr;
  hessian_[0] /= 4.0 * vars_[0] * vars_[0];
  hessian_[1] /= 4.0 * vars_[1] * vars_[0];
  hessian_[2] /= 4.0 * vars_[1] * vars_[1];
  double df = 2.0 / std::pow(rr, 4) /
              (g_e * g_e * hessian_[2] + g_g * g_g * hessian_[0] -
               2.0 * g_e * g_g * hessian_[1]);
  df *= (hessian_[0] * hessian_[2] - hessian_[1] * hessian_[1]);
  return std::array<double, 3>{beta, fstat, df};
}
}  // namespace snplib