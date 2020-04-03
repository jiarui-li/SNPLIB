#ifndef SNPLIB_UNI_MLM_H
#define SNPLIB_UNI_MLM_H

#ifndef USE_OPENBLAS
#include <mkl.h>
#else
#include <cblas.h>
#include <lapacke.h>
#endif
#include <algorithm>
#include <array>
#include <cmath>

namespace snplib {
class UniMLM {
 private:
  size_t num_samples_;
  size_t num_covariates_;
  int32_t m;
  int32_t n;
  int32_t lwork;

  const double *trait_;
  const double *lambda_;
  const double *covariates_;
  double *Q;
  double *h;
  double *p;
  double *t;
  double *r;
  double *tau;
  double *work;

  std::array<double, 2> vars_;
  std::array<double, 2> old_vars_;
  std::array<double, 2> gradients_;
  std::array<double, 2> step_;
  std::array<double, 3> hessian_;

  double f_old, f_new;

 public:
  UniMLM(const double *lambda, const double *covariates, size_t num_samples,
         size_t num_covariates);
  ~UniMLM() noexcept;
  double CalcLikelihood();
  void UpdateGradients();
  void UpdateHessian();
  double CalcLineGradient();
  void CalcInitialGuess(const double *trait);
  void CalcEMStep();
  void CalcAIStep();
  void BackupVars();
  void UpdateVars(double a);
  void CalcRes(double *res);
  void CalcBeta(double *beta);
  std::array<double, 3> CalcFTest();
  const std::array<double, 2> &GetVars() const { return vars_; }
};
}  // namespace snplib

#endif  // SNPLIB_UNI_MLM_H