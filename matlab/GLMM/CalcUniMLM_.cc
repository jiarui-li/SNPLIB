#include "../../src/genetic_variances.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto *traits = mxGetPr(prhs[0]);
  auto num_samples = static_cast<size_t>(mxGetM(prhs[0]));
  auto num_traits = static_cast<size_t>(mxGetN(prhs[0]));
  auto *covariates = mxGetPr(prhs[1]);
  auto num_covariates = static_cast<size_t>(mxGetN(prhs[1]));
  auto *lambda = mxGetPr(prhs[2]);
  auto num_threads = static_cast<size_t>(mxGetScalar(prhs[3]));
  plhs[0] = mxCreateDoubleMatrix(2, num_traits, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(num_samples, num_traits, mxREAL);
  auto *vars = mxGetPr(plhs[0]);
  auto *res = mxGetPr(plhs[1]);
  snplib::CalcUniLMM(traits, covariates, lambda, num_samples, num_covariates,
                     num_traits, vars, res, num_threads);
}
