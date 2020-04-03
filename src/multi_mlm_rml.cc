#include "multi_mlm_rml.h"

namespace snplib {
void MultiMLMRML::VectorizeTraits(const double *traits) {
  for (size_t i = 0; i < num_dims_; i++) {
    for (size_t j = 0; j < num_samples_; j++) {
      vector_traits_[j * num_dims_ + i] = traits[i * num_samples_ + j];
    }
  }
}
double MultiMLMRML::CalcLikelihood() {
  auto *v = vars_;
  LAPACKE_dtpttr(LAPACK_COL_MAJOR, 'L', p, v, Le, p);
  std::fill(Ve, Ve + num_dims_ * num_dims_, 0.0);
  for (int32_t i = 0; i < p; ++i) {
    int32_t t = p - i;
    cblas_dsyr(CblasColMajor, CblasLower, t, 1.0, v, 1, Ve + i * p + i, p);
    v += t;
  }
  LAPACKE_dtpttr(LAPACK_COL_MAJOR, 'L', p, v, Lg, p);
  std::fill(Vg, Vg + num_dims_ * num_dims_, 0.0);
  for (int32_t i = 0; i < p; ++i) {
    int32_t t = p - i;
    cblas_dsyr(CblasColMajor, CblasLower, t, 1.0, v, 1, Vg + i * p + i, p);
    v += t;
  }
  double f = 0.0;
  for (size_t l = 0; l < num_samples_; l++) {
    auto eigen = lambda_[l];
    auto *tmp_H = H + l * num_dims_ * num_dims_;
    for (size_t i = 0; i < num_dims_; i++) {
      for (size_t j = i; j < num_dims_; j++) {
        tmp_H[i * num_dims_ + j] =
            Ve[i * num_dims_ + j] + eigen * Vg[i * num_dims_ + j];
      }
    }
    LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', p, tmp_H, p);
    for (size_t i = 0; i < num_dims_; i++) {
      f += std::log(tmp_H[i * num_dims_ + i] * tmp_H[i * num_dims_ + i]);
    }
  }
  std::copy(H, H + num_samples_ * num_dims_ * num_dims_, V);
  for (size_t l = 0; l < num_samples_; l++) {
    auto *tmp_V = V + l * num_dims_ * num_dims_;
    auto *tmp_t = t + l * num_dims_;
    auto *tmp_y = vector_traits_ + l * num_dims_;
    LAPACKE_dpotri(LAPACK_COL_MAJOR, 'L', p, tmp_V, p);
    cblas_dsymv(CblasColMajor, CblasLower, p, 1.0, tmp_V, p, tmp_y, 1, 0.0,
                tmp_t, 1);
  }
  for (size_t i = 0; i < num_dims_ * num_samples_; ++i) {
    f += vector_traits_[i] * t[i];
  }
  return f / 2.0;
}
void MultiMLMRML::UpdateGradients() {
  std::fill(gradients_, gradients_ + num_dims_ * (num_dims_ + 1), 0.0);
  std::fill(Ht, Ht + num_dims_ * (num_dims_ + 1) * num_dims_ * num_samples_,
            0.0);
  for (size_t l = 0; l < num_samples_; l++) {
    auto *tmp_V = V + l * num_dims_ * num_dims_;
    auto *tmp_t = t + l * num_dims_;
    auto *tmp_ht = Ht + l * num_dims_;
    auto eigen = lambda_[l];
    for (size_t i = 0; i < num_dims_; i++) {
      for (size_t j = i + 1; j < num_dims_; j++) {
        tmp_V[j * num_dims_ + i] = tmp_V[i * num_dims_ + j];
      }
    }
    auto *gg = gradients_;
    for (size_t i = 0; i < num_dims_; ++i) {
      auto tt = num_dims_ - i;
      auto *tmp_L = Le + i * num_dims_ + i;
      for (size_t j = i; j < num_dims_; ++j) {
        auto *tmp_vv = tmp_V + j * num_dims_ + i;
        auto tmp = tmp_t[j];
        for (size_t k = 0; k < tt; k++) {
          tmp_ht[i + k] = tmp * tmp_L[k];
        }
        tmp = 0.0;
        for (size_t k = 0; k < tt; k++) {
          tmp += tmp_t[i + k] * tmp_L[k];
        }
        tmp_ht[j] += tmp;
        tmp = 0.0;
        for (size_t k = 0; k < tt; k++) {
          tmp -= tmp_t[i + k] * tmp_ht[i + k];
        }
        tmp /= 2.0;
        for (size_t k = 0; k < tt; k++) {
          tmp += tmp_L[k] * tmp_vv[k];
        }
        *gg += tmp;
        gg++;
        tmp_ht += num_dims_ * num_samples_;
      }
    }
    for (size_t i = 0; i < num_dims_; ++i) {
      auto tt = num_dims_ - i;
      auto *tmp_L = Lg + i * num_dims_ + i;
      for (size_t j = i; j < num_dims_; ++j) {
        auto *tmp_vv = tmp_V + j * num_dims_ + i;
        auto tmp = eigen * tmp_t[j];
        for (size_t k = 0; k < tt; k++) {
          tmp_ht[i + k] = tmp * tmp_L[k];
        }
        tmp = 0.0;
        for (size_t k = 0; k < tt; k++) {
          tmp += tmp_t[i + k] * tmp_L[k];
        }
        tmp_ht[j] += eigen * tmp;
        tmp = 0.0;
        for (size_t k = 0; k < tt; k++) {
          tmp -= tmp_t[i + k] * tmp_ht[i + k];
        }
        tmp /= 2.0;
        for (size_t k = 0; k < tt; k++) {
          tmp += eigen * tmp_L[k] * tmp_vv[k];
        }
        *gg += tmp;
        gg++;
        tmp_ht += num_dims_ * num_samples_;
      }
    }
  }
}
void MultiMLMRML::UpdateHessian() {
  for (size_t l = 0; l < num_samples_; l++) {
    auto *tmp_H = H + l * num_dims_ * num_dims_;
    auto *tmp_ht = Ht + l * num_dims_;
    cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans,
                CblasNonUnit, p, d, 1.0, tmp_H, p, tmp_ht, mp);
  }
  cblas_dsyrk(CblasColMajor, CblasLower, CblasTrans, d, mp, 1.0, Ht, mp, 0.0,
              hessian_, d);
}
MultiMLMRML::MultiMLMRML(const double *lambda, size_t num_samples,
                         size_t num_covariates, size_t num_dims)
    : num_samples_(num_samples),
      num_covariates_(num_covariates),
      num_para_dims_(num_dims * (num_dims + 1)),
      num_dims_(num_dims),
      lambda_(lambda) {
  l_sum = 0.0;
  for (size_t i = 0; i < num_samples_; ++i) {
    l_sum += lambda_[i];
  }
  l_sum /= num_samples_;
  m = static_cast<int32_t>(num_samples);
  p = static_cast<int32_t>(num_dims);
  mp = m * p;
  d = p * (p + 1);
  vector_traits_ = new double[num_dims * num_samples];
  Le = new double[num_dims * num_dims];
  Lg = new double[num_dims * num_dims];
  Ve = new double[num_dims * num_dims];
  Vg = new double[num_dims * num_dims];
  H = new double[num_samples * num_dims * num_dims];
  V = new double[num_samples * num_dims * num_dims];
  t = new double[num_dims * num_samples];
  Ht = new double[num_dims * (num_dims + 1) * num_dims * num_samples];
  vars_ = new double[num_para_dims_];
  old_vars_ = new double[num_para_dims_];
  gradients_ = new double[num_para_dims_];
  step_ = new double[num_para_dims_];
  hessian_ = new double[num_para_dims_ * num_para_dims_];
}
MultiMLMRML::~MultiMLMRML() noexcept {
  delete[] vector_traits_;
  delete[] Le;
  delete[] Lg;
  delete[] Ve;
  delete[] Vg;
  delete[] H;
  delete[] V;
  delete[] t;
  delete[] Ht;
  delete[] vars_;
  delete[] old_vars_;
  delete[] gradients_;
  delete[] step_;
  delete[] hessian_;
}
double MultiMLMRML::CalcLineGradient() {
  double result = 0.0;
  for (size_t i = 0; i < num_para_dims_; ++i) {
    result += gradients_[i] * step_[i];
  }
  return result;
}
void MultiMLMRML::CalcEMStep() {
  double c = 1.0 / num_samples_;
  std::fill(Le, Le + num_dims_ * num_dims_, 0.0);
  std::fill(Lg, Lg + num_dims_ * num_dims_, 0.0);
  for (size_t l = 0; l < num_samples_; l++) {
    auto *tmp_V = V + l * num_dims_ * num_dims_;
    auto *tmp_t = t + l * num_dims_;
    auto eigen = lambda_[l];
    cblas_dsyr(CblasColMajor, CblasLower, p, -1.0, tmp_t, 1, tmp_V, p);
    for (size_t i = 0; i < num_dims_; i++) {
      for (size_t j = i; j < num_dims_; j++) {
        Le[i * num_dims_ + j] += tmp_V[i * num_dims_ + j];
        Lg[i * num_dims_ + j] += eigen * tmp_V[i * num_dims_ + j];
      }
    }
  }
  for (size_t i = 0; i < num_dims_; i++) {
    for (size_t j = i; j < num_dims_; j++) {
      Le[j * num_dims_ + i] = Le[i * num_dims_ + j];
      Lg[j * num_dims_ + i] = Lg[i * num_dims_ + j];
    }
  }
  cblas_dsymm(CblasColMajor, CblasLeft, CblasLower, p, p, 1.0, Ve, p, Le, p,
              0.0, H, p);
  cblas_dsymm(CblasColMajor, CblasRight, CblasLower, p, p, c, Ve, p, H, p, 0.0,
              Le, p);
  cblas_dsymm(CblasColMajor, CblasLeft, CblasLower, p, p, 1.0, Vg, p, Lg, p,
              0.0, H, p);
  cblas_dsymm(CblasColMajor, CblasRight, CblasLower, p, p, c, Vg, p, H, p, 0.0,
              Lg, p);
  for (size_t i = 0; i < num_dims_; i++) {
    for (size_t j = i; j < num_dims_; j++) {
      Ve[i * num_dims_ + j] -= Le[i * num_dims_ + j];
      Vg[i * num_dims_ + j] -= Lg[i * num_dims_ + j];
    }
  }
  LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', p, Ve, p);
  LAPACKE_dtrttp(LAPACK_COL_MAJOR, 'L', p, Ve, p, vars_);
  LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', p, Vg, p);
  LAPACKE_dtrttp(LAPACK_COL_MAJOR, 'L', p, Vg, p, vars_ + d / 2);
}
void MultiMLMRML::CalcAIStep() {
  for (size_t i = 0; i < num_para_dims_; i++) {
    step_[i] = -gradients_[i];
  }
  LAPACKE_dposv(LAPACK_COL_MAJOR, 'L', d, 1, hessian_, d, step_, d);
}
void MultiMLMRML::CalcInitialGuess(double *res, double *vars) {
  double c = 1.0 / (num_samples_ - num_covariates_);
  cblas_dsyrk(CblasColMajor, CblasLower, CblasTrans, p, m, c, res, m, 0.0, Le,
              p);
  for (size_t i = 0; i < num_dims_; ++i) {
    for (size_t j = i; j < num_dims_; ++j) {
      auto a =
          Le[i * num_dims_ + j] / (l_sum * vars[2 * i + 1] * vars[2 * j + 1] +
                                   vars[2 * i] * vars[2 * j]);
      Ve[i * num_dims_ + j] = a * vars[2 * i] * vars[2 * j];
      Vg[i * num_dims_ + j] = a * vars[2 * i + 1] * vars[2 * j + 1];
    }
  }
  VectorizeTraits(res);
  LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', p, Ve, p);
  LAPACKE_dtrttp(LAPACK_COL_MAJOR, 'L', p, Ve, p, vars_);
  LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', p, Vg, p);
  LAPACKE_dtrttp(LAPACK_COL_MAJOR, 'L', p, Vg, p, vars_ + d / 2);
}
void MultiMLMRML::BackupVars() {
  std::copy(vars_, vars_ + num_para_dims_, old_vars_);
}
void MultiMLMRML::UpdateVars(double a) {
  for (size_t i = 0; i < num_para_dims_; ++i) {
    vars_[i] = old_vars_[i] + a * step_[i];
  }
}
void MultiMLMRML::GetSigmaE(double *sigma_e) {
  for (size_t i = 0; i < num_dims_; i++) {
    for (size_t j = i; j < num_dims_; j++) {
      sigma_e[i * num_dims_ + j] = Ve[i * num_dims_ + j];
      sigma_e[j * num_dims_ + i] = Ve[i * num_dims_ + j];
    }
  }
}
void MultiMLMRML::GetSigmaG(double *sigma_g) {
  for (size_t i = 0; i < num_dims_; i++) {
    for (size_t j = i; j < num_dims_; j++) {
      sigma_g[i * num_dims_ + j] = Vg[i * num_dims_ + j];
      sigma_g[j * num_dims_ + i] = Vg[i * num_dims_ + j];
    }
  }
}
}  // namespace snplib