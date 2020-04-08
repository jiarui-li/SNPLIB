#ifndef _SNPLIB_SRC_MULTI_LMM_H_
#define _SNPLIB_SRC_MULTI_LMM_H_

#include <algorithm>
#include <cmath>
#ifdef USE_MKL
#include <mkl.h>
#endif
#ifdef USE_OPENBLAS
#include <cblas.h>
#include <lapacke.h>
#endif

namespace snplib {
class MultiLMMREML {
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
  MultiLMMREML(const double *lambda, const double *covariates,
               size_t num_samples, size_t num_covariates, size_t num_dims);
  ~MultiLMMREML() noexcept;
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

class MultiLMMRML {
 private:
  size_t num_samples_;
  size_t num_covariates_;
  size_t num_dims_;
  size_t num_para_dims_;
  int32_t mp;
  int32_t m;
  int32_t p;
  int32_t d;
  int32_t lwork;

  double l_sum;
  const double *lambda_;
  double *vector_traits_;
  double *Le;
  double *Lg;
  double *Ve;
  double *Vg;
  double *H;
  double *V;
  double *t;
  double *Ht;

  double *vars_;
  double *old_vars_;
  double *gradients_;
  double *step_;
  double *hessian_;
  void VectorizeTraits(const double *traits);

 public:
  MultiLMMRML(const double *lambda, size_t num_samples, size_t num_covariates_,
              size_t num_dims);
  ~MultiLMMRML() noexcept;
  double CalcLikelihood();
  void UpdateGradients();
  void UpdateHessian();
  double CalcLineGradient();
  void CalcInitialGuess(double *res, double *vars);
  void CalcEMStep();
  void CalcAIStep();
  void BackupVars();
  void UpdateVars(double a);
  void GetSigmaE(double *sigma_e);
  void GetSigmaG(double *sigma_g);
};
}  // namespace snplib

#endif  //_SNPLIB_SRC_MULTI_LMM_H_
