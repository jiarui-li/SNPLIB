#include "../../src/genetic_variances.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto *traits = mxGetPr(prhs[0]);
  auto num_samples = static_cast<size_t>(mxGetM(prhs[0]));
  auto num_traits = static_cast<size_t>(mxGetN(prhs[0]));
  auto *covariates = mxGetPr(prhs[1]);
  auto num_covariates = static_cast<size_t>(mxGetN(prhs[1]));
  auto *lambda = mxGetPr(prhs[2]);
  auto *res = mxGetPr(prhs[3]);
  auto *vars = mxGetPr(prhs[4]);
  auto num_dims = static_cast<size_t>(mxGetScalar(prhs[5]));
  auto num_threads = static_cast<size_t>(mxGetScalar(prhs[6]));
  num_traits /= num_dims;
  mwSize dims[3] = {num_dims, num_dims, num_traits};
  plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  auto *sigma_e = mxGetPr(plhs[0]);
  auto *sigma_g = mxGetPr(plhs[1]);
  snplib::CalcMultiLMM_REML(covariates, lambda, traits, res, vars, num_samples,
                            num_covariates, num_dims, num_traits, sigma_e,
                            sigma_g, num_threads);
}