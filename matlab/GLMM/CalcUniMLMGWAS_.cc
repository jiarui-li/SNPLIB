#include "../../src/gwas.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto *trait = mxGetPr(prhs[0]);
  auto *geno = reinterpret_cast<uint8_t *>(mxGetData(prhs[1]));
  auto num_snps = static_cast<size_t>(mxGetN(prhs[1]));
  auto *covariates = mxGetPr(prhs[2]);
  auto num_samples = static_cast<size_t>(mxGetM(prhs[2]));
  auto num_covariates = static_cast<size_t>(mxGetN(prhs[2]));
  auto *lambda = mxGetPr(prhs[3]);
  auto num_threads = static_cast<size_t>(mxGetScalar(prhs[4]));
  plhs[0] = mxCreateDoubleMatrix(1, num_snps, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1, num_snps, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(1, num_snps, mxREAL);
  auto *betas = mxGetPr(plhs[0]);
  auto *fstats = mxGetPr(plhs[1]);
  auto *dfs = mxGetPr(plhs[2]);
  snplib::CalcUniLMMGWAS(trait, geno, covariates, lambda, num_samples,
                         num_covariates, num_snps, betas, fstats, dfs,
                         num_threads);
}