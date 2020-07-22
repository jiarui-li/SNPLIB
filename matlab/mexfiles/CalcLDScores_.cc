#include "../../src/statistics.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto *geno = reinterpret_cast<uint8_t *>(mxGetData(prhs[0]));
  auto num_snps = static_cast<size_t>(mxGetN(prhs[0]));
  auto *bp = reinterpret_cast<int32_t *>(mxGetData(prhs[1]));
  auto *af = mxGetPr(prhs[2]);
  auto num_samples = static_cast<size_t>(mxGetScalar(prhs[3]));
  auto window_size = static_cast<size_t>(mxGetScalar(prhs[4]));
  auto r2_threshold = mxGetScalar(prhs[5]);
  auto num_threads = static_cast<size_t>(mxGetScalar(prhs[6]));
  plhs[0] = mxCreateDoubleMatrix(num_snps, 1, mxREAL);
  auto *ldcv = mxGetPr(plhs[0]);
  snplib::CalcLDscores(geno, bp, af, num_snps, num_samples, window_size,
                       r2_threshold, ldcv, num_threads);
}