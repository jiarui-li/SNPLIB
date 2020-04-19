#include "../../src/statistics.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto *geno = (uint8_t *)mxGetData(prhs[0]);
  auto num_snps = static_cast<size_t>(mxGetN(prhs[0]));
  auto *covariates = mxGetPr(prhs[1]);
  auto num_samples = static_cast<size_t>(mxGetM(prhs[1]));
  auto num_covariates = static_cast<size_t>(mxGetN(prhs[1]));
  auto num_threads = static_cast<size_t>(mxGetScalar(prhs[2]));
  plhs[0] = mxCreateDoubleMatrix(num_snps, 1, mxREAL);
  auto *min_maf = mxGetPr(plhs[0]);
  snplib::CalcAdjustedMAF(geno, num_samples, num_snps, covariates,
                          num_covariates, min_maf, num_threads);
}
