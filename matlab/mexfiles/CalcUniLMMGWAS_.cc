#include "../../src/gwas.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto *geno = reinterpret_cast<uint8_t *>(mxGetData(prhs[0]));
  auto num_snps = static_cast<size_t>(mxGetN(prhs[0]));
  auto *covariates = mxGetPr(prhs[1]);
  auto num_samples = static_cast<size_t>(mxGetM(prhs[1]));
  auto num_covariates = static_cast<size_t>(mxGetN(prhs[1]));
  auto *trait = mxGetPr(prhs[2]);
  auto *lambda = mxGetPr(prhs[3]);
  auto *V = mxGetPr(prhs[4]);
  auto num_threads = static_cast<size_t>(mxGetScalar(prhs[5]));
  plhs[0] = mxCreateDoubleMatrix(num_snps, 1, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(num_snps, 1, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(num_snps, 1, mxREAL);
  auto *betas = mxGetPr(plhs[0]);
  auto *fstats = mxGetPr(plhs[1]);
  auto *dfs = mxGetPr(plhs[2]);
  snplib::CalcUniLMMGWAS(geno, num_samples, num_snps, lambda, V, covariates,
                         num_covariates, trait, betas, fstats, dfs,
                         num_threads);
}