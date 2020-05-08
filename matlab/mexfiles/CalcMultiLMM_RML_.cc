#include "../../src/genetic_variances.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto *lambda = mxGetPr(prhs[0]);
  auto *res = mxGetPr(prhs[1]);
  auto num_samples = static_cast<size_t>(mxGetM(prhs[1]));
  auto *vars = mxGetPr(prhs[2]);
  auto num_traits = static_cast<size_t>(mxGetN(prhs[2]));
  auto num_covariates = static_cast<size_t>(mxGetScalar(prhs[3]));
  auto num_dims = static_cast<size_t>(mxGetScalar(prhs[4]));
  auto num_threads = static_cast<size_t>(mxGetScalar(prhs[5]));
  num_traits /= num_dims;
  mwSize dims[3] = {num_dims, num_dims, num_traits};
  plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  auto *sigma_e = mxGetPr(plhs[0]);
  auto *sigma_g = mxGetPr(plhs[1]);
  snplib::CalcMultiLMM_RML(lambda, res, vars, num_samples, num_covariates,
                           num_dims, num_traits, sigma_e, sigma_g, num_threads);
}