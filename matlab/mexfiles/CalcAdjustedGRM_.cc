#include "../../src/adjusted_grm.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto *geno = reinterpret_cast<uint8_t *>(mxGetData(prhs[0]));
  auto num_snps = static_cast<size_t>(mxGetN(prhs[0]));
  auto *covariates = mxGetPr(prhs[1]);
  auto num_samples = static_cast<size_t>(mxGetM(prhs[1]));
  auto num_covariates = static_cast<size_t>(mxGetN(prhs[1]));
  auto num_threads = static_cast<size_t>(mxGetScalar(prhs[2]));
  plhs[0] = mxCreateDoubleMatrix(num_samples, num_samples, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(num_samples, 1, mxREAL);
  auto *matrix = mxGetPr(plhs[0]);
  auto *gcta_diag = mxGetPr(plhs[1]);
  snplib::CalcAdjustedGRM(geno, num_samples, num_snps, covariates,
                          num_covariates, matrix, gcta_diag, num_threads);
}