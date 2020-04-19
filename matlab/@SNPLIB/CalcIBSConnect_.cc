#include "../../src/relationships.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto *src_geno = (uint8_t *)mxGetData(prhs[0]);
  auto num_snps = static_cast<size_t>(mxGetN(prhs[0]));
  auto num_src_samples = static_cast<size_t>(mxGetScalar(prhs[1]));
  auto *dest_geno = (uint8_t *)mxGetData(prhs[2]);
  auto num_dest_samples = static_cast<size_t>(mxGetScalar(prhs[3]));
  auto num_threads = static_cast<size_t>(mxGetScalar(prhs[4]));
  plhs[0] = mxCreateDoubleMatrix(num_dest_samples, 1, mxREAL);
  auto *connect = mxGetPr(plhs[0]);
  snplib::CalcIBSConnection(src_geno, num_src_samples, dest_geno,
                            num_dest_samples, num_snps, connect, num_threads);
}
