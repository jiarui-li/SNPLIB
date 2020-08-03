#include <algorithm>
#include <cmath>
#include <cstdint>

#include "blas.h"
#include "lapacke.h"
#include "mex.h"

namespace {
template <class T>
class LineSearch {
 private:
  T &worker_;
  double f_new;
  double f_old, g_old;
  double a_h, a_l;
  double f_h, f_l;
  double g_l;

  double Zoom() {
    double a_j;
    for (size_t i = 0; i < 10; ++i) {
      auto t = a_h - a_l;
      auto a = (f_h - f_l) - t * g_l;
      a /= t * t;
      a_j = a_l - g_l / a / 2.0;
      worker_.UpdateVars(a_j);
      double f_j = worker_.CalcLikelihood();
      if (a_j == a_l || i == 9) {
        worker_.UpdateGradients();
        f_new = f_j;
        return a_j;
      }
      if ((f_j - f_old) > 1e-4 * a_j * g_old || f_j >= f_l) {
        a_h = a_j;
        f_h = f_j;
      } else {
        worker_.UpdateGradients();
        double g_j = worker_.CalcLineGradient();
        if (std::fabs(g_j) <= -0.9 * g_old) {
          f_new = f_j;
          return a_j;
        }
        if (g_j * (a_h - a_l) >= 0) {
          a_h = a_l;
          f_h = f_l;
        }
        a_l = a_j;
        f_l = f_j;
        g_l = g_j;
      }
    }
    return a_j;
  }

 public:
  LineSearch(T &worker) : worker_(worker) {}
  ~LineSearch() {}
  double Search(double f) {
    g_old = worker_.CalcLineGradient();
    f_old = f;
    double a_1 = 1.0;
    double a_0 = 0.0;
    double f_0 = f_old;
    double g_0 = g_old;
    size_t num_iter = 0;
    while (true) {
      worker_.UpdateVars(a_1);
      double f_1 = worker_.CalcLikelihood();
      if ((f_1 - f_old) > 1e-4 * a_1 * g_old || (num_iter > 0 && f_1 >= f_0)) {
        a_h = a_1;
        f_h = f_1;
        a_l = a_0;
        f_l = f_0;
        g_l = g_0;
        return Zoom();
      }
      worker_.UpdateGradients();
      double g_1 = worker_.CalcLineGradient();
      if (std::fabs(g_1) <= -0.9 * g_old) {
        f_new = f_1;
        return a_1;
      }
      if (g_1 >= 0) {
        a_h = a_0;
        a_l = a_1;
        f_h = f_0;
        f_l = f_1;
        g_l = g_1;
        return Zoom();
      }
      a_0 = a_1;
      f_0 = f_1;
      g_0 = g_1;
      a_1 *= 1.618;
      num_iter++;
    }
  }
  double GetFNew() const { return f_new; }
};
void VectorizeTraits(const double *traits, size_t num_dims, size_t num_samples,
                     double *vector_traits) {
  for (size_t i = 0; i < num_dims; i++) {
    for (size_t j = 0; j < num_samples; j++) {
      vector_traits[j * num_dims + i] = traits[i * num_samples + j];
    }
  }
}
class MultiLMM_REML {
 private:
  size_t num_samples_;
  size_t num_covariates_;
  size_t num_dims_;
  size_t num_para_dims_;
  ptrdiff_t m;
  ptrdiff_t mp;
  ptrdiff_t np;
  ptrdiff_t p;
  ptrdiff_t d;
  ptrdiff_t lwork;

  double l_sum;
  const double *lambda_;
  const double *covariates_;
  double *vector_traits_;
  double *Le;
  double *Lg;
  double *Ve;
  double *Vg;
  double *Q;
  double *H;
  double *P;
  double *t;
  double *r;
  double *Ht;
  double *QHt;
  double *tau;
  double *work;

  double *vars_;
  double *old_vars_;
  double *gradients_;
  double *step_;
  double *hessian_;
  void FillQ() {
    std::fill(Q, Q + num_dims_ * num_dims_ * num_covariates_ * num_samples_,
              0.0);
    for (size_t i = 0; i < num_covariates_; ++i) {
      auto *tmp_c = covariates_ + i * num_samples_;
      for (size_t j = 0; j < num_samples_; ++j) {
        auto tmp = tmp_c[j];
        for (size_t k = 0; k < num_dims_; ++k) {
          Q[(i * num_dims_ + k) * num_dims_ * num_samples_ +
            (j * num_dims_ + k)] = tmp;
        }
      }
    }
  }

 public:
  MultiLMM_REML(const double *lambda, const double *covariates,
                size_t num_samples, size_t num_covariates, size_t num_dims)
      : num_samples_(num_samples),
        num_covariates_(num_covariates),
        num_dims_(num_dims),
        num_para_dims_(num_dims * (num_dims + 1)),
        lambda_(lambda),
        covariates_(covariates) {
    l_sum = 0.0;
    for (size_t i = 0; i < num_samples_; ++i) {
      l_sum += lambda_[i];
    }
    l_sum /= num_samples_;
    m = static_cast<ptrdiff_t>(num_samples);
    mp = static_cast<ptrdiff_t>(num_dims * num_samples);
    np = static_cast<ptrdiff_t>(num_dims * num_covariates);
    p = static_cast<ptrdiff_t>(num_dims);
    d = p * (p + 1);
    vector_traits_ = new double[num_dims * num_samples];
    Le = new double[num_dims * num_dims];
    Lg = new double[num_dims * num_dims];
    Ve = new double[num_dims * num_dims];
    Vg = new double[num_dims * num_dims];
    Q = new double[num_dims * num_samples * num_dims * num_covariates];
    H = new double[num_samples * num_dims * num_dims];
    P = new double[num_samples * num_dims * num_dims];
    t = new double[num_dims * num_samples];
    r = new double[num_dims * num_covariates];
    Ht = new double[num_dims * (num_dims + 1) * num_dims * num_samples];
    QHt = new double[num_dims * (num_dims + 1) * num_dims * num_covariates];
    tau = new double[num_dims * num_covariates];
    double w1;
    LAPACKE_dgeqrf_work(LAPACK_COL_MAJOR, mp, np, nullptr, mp, nullptr, &w1,
                        -1);
    double w2;
    LAPACKE_dorgqr_work(LAPACK_COL_MAJOR, mp, np, np, nullptr, mp, nullptr, &w2,
                        -1);
    lwork = static_cast<ptrdiff_t>(w1 > w2 ? w1 : w2);
    work = new double[lwork];
    vars_ = new double[num_para_dims_];
    old_vars_ = new double[num_para_dims_];
    gradients_ = new double[num_para_dims_];
    step_ = new double[num_para_dims_];
    hessian_ = new double[num_para_dims_ * num_para_dims_];
  }
  ~MultiLMM_REML() {
    delete[] vector_traits_;
    delete[] Le;
    delete[] Lg;
    delete[] Ve;
    delete[] Vg;
    delete[] Q;
    delete[] H;
    delete[] P;
    delete[] t;
    delete[] r;
    delete[] Ht;
    delete[] QHt;
    delete[] tau;
    delete[] work;
    delete[] vars_;
    delete[] old_vars_;
    delete[] gradients_;
    delete[] step_;
    delete[] hessian_;
  }
  double CalcLikelihood() {
    auto *v = vars_;
    LAPACKE_dtpttr(LAPACK_COL_MAJOR, 'L', p, v, Le, p);
    std::fill(Ve, Ve + num_dims_ * num_dims_, 0.0);
    double alpha = 1.0;
    ptrdiff_t inc = 1;
    for (ptrdiff_t i = 0; i < p; ++i) {
      ptrdiff_t t = p - i;
      dsyr("L", &t, &alpha, v, &inc, Ve + i * p + i, &p);
      v += t;
    }
    LAPACKE_dtpttr(LAPACK_COL_MAJOR, 'L', p, v, Lg, p);
    std::fill(Vg, Vg + num_dims_ * num_dims_, 0.0);
    for (ptrdiff_t i = 0; i < p; ++i) {
      ptrdiff_t t = p - i;
      dsyr("L", &t, &alpha, v, &inc, Vg + i * p + i, &p);
      v += t;
    }
    double f = 0.0;
    FillQ();
    for (size_t l = 0; l < num_samples_; l++) {
      auto eigen = lambda_[l];
      auto *tmp_H = H + l * num_dims_ * num_dims_;
      auto *tmp_Q = Q + l * num_dims_;
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
      dtrsm("L", "L", "N", "N", &p, &np, &alpha, tmp_H, &p, tmp_Q, &mp);
    }
    LAPACKE_dgeqrf_work(LAPACK_COL_MAJOR, mp, np, Q, mp, tau, work, lwork);
    for (size_t i = 0; i < num_dims_ * num_covariates_; ++i) {
      f += std::log(Q[i * num_dims_ * num_samples_ + i] *
                    Q[i * num_dims_ * num_samples_ + i]);
    }
    LAPACKE_dorgqr_work(LAPACK_COL_MAJOR, mp, np, np, Q, mp, tau, work, lwork);
    std::copy(H, H + num_dims_ * num_dims_ * num_samples_, P);
    double beta = 0.0;
    for (size_t l = 0; l < num_samples_; l++) {
      auto *tmp_H = H + l * num_dims_ * num_dims_;
      auto *tmp_P = P + l * num_dims_ * num_dims_;
      auto *tmp_Q = Q + l * num_dims_;
      auto *tmp_y = vector_traits_ + l * num_dims_;
      auto *tmp_t = t + l * num_dims_;
      LAPACKE_dpotri(LAPACK_COL_MAJOR, 'L', p, tmp_P, p);
      dtrsm("L", "L", "T", "N", &p, &np, &alpha, tmp_H, &p, tmp_Q, &mp);
      dsymv("L", &p, &alpha, tmp_P, &p, tmp_y, &inc, &beta, tmp_t, &inc);
    }
    dgemv("T", &mp, &np, &alpha, Q, &mp, vector_traits_, &inc, &beta, r, &inc);
    alpha = -1.0;
    beta = 1.0;
    dgemv("N", &mp, &np, &alpha, Q, &mp, r, &inc, &beta, t, &inc);
    for (size_t i = 0; i < num_dims_ * num_samples_; ++i) {
      f += t[i] * vector_traits_[i];
    }
    return f / 2.0;
  }
  void UpdateGradients() {
    std::fill(gradients_, gradients_ + num_dims_ * (num_dims_ + 1), 0.0);
    std::fill(Ht, Ht + num_dims_ * (num_dims_ + 1) * num_dims_ * num_samples_,
              0.0);
    double alpha = -1.0;
    double beta = 1.0;
    for (size_t l = 0; l < num_samples_; l++) {
      auto *tmp_P = P + l * num_dims_ * num_dims_;
      auto *tmp_Q = Q + l * num_dims_;
      auto *tmp_t = t + l * num_dims_;
      auto *tmp_ht = Ht + l * num_dims_;
      auto eigen = lambda_[l];
      dsyrk("L", "N", &p, &np, &alpha, tmp_Q, &mp, &beta, tmp_P, &p);
      for (size_t i = 0; i < num_dims_; ++i) {
        for (size_t j = i + 1; j < num_dims_; ++j) {
          tmp_P[j * num_dims_ + i] = tmp_P[i * num_dims_ + j];
        }
      }
      auto *gg = gradients_;
      for (size_t i = 0; i < num_dims_; ++i) {
        auto tt = num_dims_ - i;
        auto *tmp_L = Le + i * num_dims_ + i;
        for (size_t j = i; j < num_dims_; ++j) {
          auto *tmp_p = tmp_P + j * num_dims_ + i;
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
            tmp += tmp_L[k] * tmp_p[k];
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
          auto *tmp_p = tmp_P + j * num_dims_ + i;
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
            tmp += eigen * tmp_L[k] * tmp_p[k];
          }
          *gg += tmp;
          gg++;
          tmp_ht += num_dims_ * num_samples_;
        }
      }
    }
  }
  void UpdateHessian() {
    double alpha = 1.0;
    double beta = 0.0;
    dgemm("T", "N", &np, &d, &mp, &alpha, Q, &mp, Ht, &mp, &beta, QHt, &np);
    for (size_t l = 0; l < num_samples_; l++) {
      auto *tmp_H = H + l * num_dims_ * num_dims_;
      auto *tmp_ht = Ht + l * num_dims_;
      dtrsm("L", "L", "N", "N", &p, &d, &alpha, tmp_H, &p, tmp_ht, &mp);
    }
    dsyrk("L", "T", &d, &mp, &alpha, Ht, &mp, &beta, hessian_, &d);
    alpha = -1.0;
    beta = 1.0;
    dsyrk("L", "T", &d, &np, &alpha, QHt, &np, &beta, hessian_, &d);
  }
  double CalcLineGradient() {
    double result = 0.0;
    for (size_t i = 0; i < num_para_dims_; ++i) {
      result += gradients_[i] * step_[i];
    }
    return result;
  }
  void CalcInitialGuess(const double *traits, const double *res,
                        const double *vars) {
    double c = 1.0 / (num_samples_ - num_covariates_);
    double beta = 0.0;
    dsyrk("L", "T", &p, &m, &c, res, &m, &beta, Le, &p);
    for (size_t i = 0; i < num_dims_; ++i) {
      for (size_t j = i; j < num_dims_; ++j) {
        auto a =
            Le[i * num_dims_ + j] / (l_sum * vars[2 * i + 1] * vars[2 * j + 1] +
                                     vars[2 * i] * vars[2 * j]);
        Ve[i * num_dims_ + j] = a * vars[2 * i] * vars[2 * j];
        Vg[i * num_dims_ + j] = a * vars[2 * i + 1] * vars[2 * j + 1];
      }
    }
    VectorizeTraits(traits, num_dims_, num_samples_, vector_traits_);
    LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', p, Ve, p);
    LAPACKE_dtrttp(LAPACK_COL_MAJOR, 'L', p, Ve, p, vars_);
    LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', p, Vg, p);
    LAPACKE_dtrttp(LAPACK_COL_MAJOR, 'L', p, Vg, p, vars_ + d / 2);
  }
  void CalcEMStep() {
    double c = 1.0 / num_samples_;
    std::fill(Le, Le + num_dims_ * num_dims_, 0.0);
    std::fill(Lg, Lg + num_dims_ * num_dims_, 0.0);
    double alpha = -1.0;
    double beta = 1.0;
    ptrdiff_t inc = 1;
    for (size_t l = 0; l < num_samples_; l++) {
      auto *tmp_P = P + l * num_dims_ * num_dims_;
      auto *tmp_Q = Q + l * num_dims_;
      auto *tmp_t = t + l * num_dims_;
      auto eigen = lambda_[l];
      dsyrk("L", "N", &p, &np, &alpha, tmp_Q, &mp, &beta, tmp_P, &p);
      dsyr("L", &p, &alpha, tmp_t, &inc, tmp_P, &p);
      for (size_t i = 0; i < num_dims_; i++) {
        for (size_t j = i; j < num_dims_; j++) {
          Le[i * num_dims_ + j] += tmp_P[i * num_dims_ + j];
          Lg[i * num_dims_ + j] += eigen * tmp_P[i * num_dims_ + j];
        }
      }
    }
    for (size_t i = 0; i < num_dims_; i++) {
      for (size_t j = i; j < num_dims_; j++) {
        Le[j * num_dims_ + i] = Le[i * num_dims_ + j];
        Lg[j * num_dims_ + i] = Lg[i * num_dims_ + j];
      }
    }
    alpha = 1.0;
    beta = 0.0;
    dsymm("L", "L", &p, &p, &alpha, Ve, &p, Le, &p, &beta, H, &p);
    dsymm("R", "L", &p, &p, &c, Ve, &p, H, &p, &beta, Le, &p);
    dsymm("L", "L", &p, &p, &alpha, Vg, &p, Lg, &p, &beta, H, &p);
    dsymm("R", "L", &p, &p, &c, Vg, &p, H, &p, &beta, Lg, &p);
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
  void CalcAIStep() {
    for (size_t i = 0; i < num_para_dims_; i++) {
      step_[i] = -gradients_[i];
    }
    LAPACKE_dposv(LAPACK_COL_MAJOR, 'L', d, 1, hessian_, d, step_, d);
  }
  void BackupVars() { std::copy(vars_, vars_ + num_para_dims_, old_vars_); }
  void UpdateVars(double a) {
    for (size_t i = 0; i < num_para_dims_; ++i) {
      vars_[i] = old_vars_[i] + a * step_[i];
    }
  }
  void GetSigmaE(double *sigma_e) const {
    for (size_t i = 0; i < num_dims_; i++) {
      for (size_t j = i; j < num_dims_; j++) {
        sigma_e[i * num_dims_ + j] = Ve[i * num_dims_ + j];
        sigma_e[j * num_dims_ + i] = Ve[i * num_dims_ + j];
      }
    }
  }

  void GetSigmaG(double *sigma_g) const {
    for (size_t i = 0; i < num_dims_; i++) {
      for (size_t j = i; j < num_dims_; j++) {
        sigma_g[i * num_dims_ + j] = Vg[i * num_dims_ + j];
        sigma_g[j * num_dims_ + i] = Vg[i * num_dims_ + j];
      }
    }
  }
};
}  // namespace

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto *traits = mxGetPr(prhs[0]);
  auto num_samples = static_cast<size_t>(mxGetM(prhs[0]));
  auto num_dims = static_cast<size_t>(mxGetN(prhs[0]));
  auto *covariates = mxGetPr(prhs[1]);
  auto num_covariates = static_cast<size_t>(mxGetN(prhs[1]));
  auto *lambda = mxGetPr(prhs[2]);
  auto *res = mxGetPr(prhs[3]);
  auto *vars = mxGetPr(prhs[4]);
  plhs[0] = mxCreateDoubleMatrix(num_dims, num_dims, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(num_dims, num_dims, mxREAL);
  auto *sigma_e = mxGetPr(plhs[0]);
  auto *sigma_g = mxGetPr(plhs[1]);
  MultiLMM_REML worker(lambda, covariates, num_samples, num_covariates,
                       num_dims);
  LineSearch<MultiLMM_REML> searcher(worker);
  mexPrintf("Calculating Initial Guess...\n");
  worker.CalcInitialGuess(traits, res, vars);
  double a;
  double f_old = worker.CalcLikelihood();
  mexPrintf("The likelihood for the initial guess is: %f \n", f_old);
  mexPrintf("Starting EM iterations...\n");
  for (size_t i = 0; i < 15; ++i) {
    worker.CalcEMStep();
    f_old = worker.CalcLikelihood();
    mexPrintf("iteration-%d: log-likelihood %f \n", i, f_old);
  }
  worker.UpdateGradients();
  double f_new = f_old;
  mexPrintf("Starting AI iterations...\n");
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
    mexPrintf("iteration-%d: log-likelihood %f \n", l, f_new);
  }
  mexPrintf("Finished!\n");
  worker.GetSigmaE(sigma_e);
  worker.GetSigmaG(sigma_g);
}