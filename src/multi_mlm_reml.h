#ifndef SNPLIB_MULTI_MLM_REML_H
#define SNPLIB_MULTI_MLM_REML_H

#ifndef USE_OPENBLAS
#include <mkl.h>
#else
#include <cblas.h>
#include <lapacke.h>
#endif
#include <algorithm>
#include <cmath>

namespace snplib {
class MultiMLMREML {
 private:
  size_t num_samples_;
  size_t num_covariates_;
  size_t num_dims_;
  size_t num_para_dims_;
  int32_t m;
  int32_t mp;
  int32_t np;
  int32_t p;
  int32_t d;
  int32_t lwork;

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
  void FillQ();
  void VectorizeTraits(const double *traits);

 public:
  MultiMLMREML(const double *lambda, const double *covariates,
               size_t num_samples, size_t num_covariates, size_t num_dims);
  ~MultiMLMREML() noexcept;
  double CalcLikelihood();
  void UpdateGradients();
  void UpdateHessian();
  double CalcLineGradient();
  void CalcInitialGuess(double *traits, double *res, double *vars);
  void CalcEMStep();
  void CalcAIStep();
  void BackupVars();
  void UpdateVars(double a);
  void CalcRes(double *Res);
  void GetSigmaE(double *sigma_e);
  void GetSigmaG(double *sigma_g);
};
}  // namespace snplib

#endif  // SNPLIB_MULTI_MLM_REML_H
