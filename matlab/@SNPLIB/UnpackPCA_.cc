#include "../../src/data_manage.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto *geno = (uint8_t *)mxGetData(prhs[0]);
  auto num_snps = static_cast<size_t>(mxGetN(prhs[0]));
  auto *af = mxGetPr(prhs[1]);
  auto num_samples = static_cast<size_t>(mxGetScalar(prhs[2]));
  plhs[0] = mxCreateDoubleMatrix(num_samples, num_snps, mxREAL);
  auto *geno_d = mxGetPr(plhs[0]);
  snplib::UnpackGRMGeno(geno, af, num_samples, num_snps, geno_d);
}
