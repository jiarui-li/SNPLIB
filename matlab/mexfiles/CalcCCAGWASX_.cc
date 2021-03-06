#include "../../src/gwas.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto *geno = reinterpret_cast<uint8_t *>(mxGetData(prhs[0]));
  auto num_snps = static_cast<size_t>(mxGetN(prhs[0]));
  auto *trait = mxGetPr(prhs[1]);
  auto num_samples = static_cast<size_t>(mxGetM(prhs[1]));
  auto num_dims = static_cast<size_t>(mxGetN(prhs[1]));
  auto *sex = mxGetPr(prhs[2]);
  auto num_threads = static_cast<size_t>(mxGetScalar(prhs[3]));
  plhs[0] = mxCreateDoubleMatrix(num_dims, num_snps, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(num_snps, 1, mxREAL);
  auto *betas = mxGetPr(plhs[0]);
  auto *rho2 = mxGetPr(plhs[1]);
  snplib::CalcCCAGWASX(geno, num_samples, num_snps, trait, sex, num_dims, betas,
                       rho2, num_threads);
}