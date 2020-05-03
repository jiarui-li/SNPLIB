#ifndef SNPLIB_SRC_MULTI_LMM_H_
#define SNPLIB_SRC_MULTI_LMM_H_

#include <algorithm>
#include <array>
#include <cmath>

#include "math_lib.h"

namespace snplib {
class MultiLMM_REML {
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

 public:
  MultiLMM_REML(const double *lambda, const double *covariates,
                size_t num_samples, size_t num_covariates, size_t num_dims);
  ~MultiLMM_REML() noexcept;
  double CalcLikelihood();
  void UpdateGradients();
  void UpdateHessian();
  double CalcLineGradient();
  void CalcInitialGuess(const double *traits, const double *res,
                        const double *vars);
  void CalcEMStep();
  void CalcAIStep();
  void BackupVars();
  void UpdateVars(double a);
  void CalcRes(double *Res);
  void GetSigmaE(double *sigma_e);
  void GetSigmaG(double *sigma_g);
};
class MultiLMM_RML {
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

 public:
  MultiLMM_RML(const double *lambda, size_t num_samples, size_t num_covariates,
               size_t num_dims);
  ~MultiLMM_RML() noexcept;
  double CalcLikelihood();
  void UpdateGradients();
  void UpdateHessian();
  double CalcLineGradient();
  void CalcInitialGuess(const double *res, const double *vars);
  void CalcEMStep();
  void CalcAIStep();
  void BackupVars();
  void UpdateVars(double a);
  void GetSigmaE(double *sigma_e);
  void GetSigmaG(double *sigma_g);
};
}  // namespace snplib

#endif  // SNPLIB_SRC_MULTI_LMM_H_
