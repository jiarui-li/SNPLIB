#include "../../src/statistics.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto *geno = (uint8_t *)mxGetData(prhs[0]);
  size_t num_snps = static_cast<size_t>(mxGetN(prhs[0]));
  size_t num_samples = static_cast<size_t>(mxGetScalar(prhs[1]));
  plhs[0] = mxCreateDoubleMatrix(num_snps, 1, mxREAL);
  auto *af = mxGetPr(plhs[0]);
  snplib::CalcAlleleFrequencies(geno, num_samples, num_snps, af);
}