#include "../../src/gwas.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto *geno = reinterpret_cast<uint8_t *>(mxGetData(prhs[0]));
  auto num_snps = static_cast<size_t>(mxGetN(prhs[0]));
  auto *scores = mxGetPr(prhs[1]);
  auto num_samples = static_cast<size_t>(mxGetM(prhs[1]));
  auto num_dims = static_cast<size_t>(mxGetN(prhs[1]));
  auto *betas = mxGetPr(prhs[2]);
  auto *sex = mxGetPr(prhs[3]);
  auto num_threads = static_cast<size_t>(mxGetScalar(prhs[4]));
  plhs[0] = mxCreateDoubleMatrix(num_snps, 1, mxREAL);
  auto *rho = mxGetPr(plhs[0]);
  snplib::CalcCCAReplicationX(geno, num_samples, num_snps, scores, betas, sex,
                              num_dims, rho, num_threads);
}