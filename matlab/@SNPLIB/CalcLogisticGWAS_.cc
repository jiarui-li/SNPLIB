#include "../../src/gwas.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto *geno = (uint8_t *)mxGetData(prhs[0]);
  auto num_snps = static_cast<size_t>(mxGetN(prhs[0]));
  auto *trait = mxGetPr(prhs[1]);
  auto num_samples = static_cast<size_t>(mxGetM(prhs[1]));
  auto *covariates = mxGetPr(prhs[2]);
  auto num_covariates = static_cast<size_t>(mxGetN(prhs[2]));
  auto num_threads = static_cast<size_t>(mxGetScalar(prhs[3]));
  plhs[0] = mxCreateDoubleMatrix(1, num_snps, mxREAL);
  auto *chi2stat = mxGetPr(plhs[0]);
  snplib::CalcLogisticGWAS(trait, covariates, geno, num_samples, num_covariates,
                           num_snps, chi2stat, num_threads);
}