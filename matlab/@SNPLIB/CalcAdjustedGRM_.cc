#include "../../src/relationships.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto *geno = (uint8_t *)mxGetData(prhs[0]);
  size_t num_snps = static_cast<size_t>(mxGetN(prhs[0]));
  auto *covariates = mxGetPr(prhs[1]);
  auto num_covariates = mxGetN(prhs[1]);
  auto num_samples = mxGetM(prhs[1]);
  auto num_threads = static_cast<size_t>(mxGetScalar(prhs[2]));
  plhs[0] = mxCreateDoubleMatrix(num_samples, num_samples, mxREAL);
  auto *matrix = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(num_samples, 1, mxREAL);
  auto *gcta_diag = mxGetPr(plhs[1]);
  snplib::CalcAdjustedGRM(geno, num_samples, num_snps, covariates,
                          num_covariates, matrix, gcta_diag, num_threads);
}