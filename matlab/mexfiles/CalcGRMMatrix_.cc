#include "../../src/grm.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto *geno = reinterpret_cast<uint8_t *>(mxGetData(prhs[0]));
  auto num_snps = static_cast<size_t>(mxGetN(prhs[0]));
  auto *af = mxGetPr(prhs[1]);
  auto num_samples = static_cast<size_t>(mxGetScalar(prhs[2]));
  auto num_threads = static_cast<size_t>(mxGetScalar(prhs[3]));
  plhs[0] = mxCreateDoubleMatrix(num_samples, num_samples, mxREAL);
  auto *matrix = mxGetPr(plhs[0]);
  snplib::CalcGRMMatrix(geno, af, num_samples, num_snps, matrix, num_threads);
}