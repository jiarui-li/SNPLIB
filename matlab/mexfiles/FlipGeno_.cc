#include "../../src/data_manage.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto *geno = reinterpret_cast<uint8_t *>(mxGetData(prhs[0]));
  auto num_samples = static_cast<size_t>(mxGetScalar(prhs[1]));
  auto *index = reinterpret_cast<int32_t *>(mxGetData(prhs[2]));
  auto num_snps = static_cast<size_t>(mxGetM(prhs[2]));
  snplib::FlipGeno(geno, num_samples, num_snps, index);
}