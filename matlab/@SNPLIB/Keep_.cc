#include "../../src/base.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto *src_geno = (uint8_t *)mxGetData(prhs[0]);
  auto num_snps = static_cast<size_t>(mxGetN(prhs[0]));
  auto num_src_samples = static_cast<size_t>(mxGetScalar(prhs[1]));
  auto num_dest_samples = static_cast<size_t>(mxGetM(prhs[2]));
  auto *sample_idx = (int32_t *)mxGetData(prhs[2]);
  auto num_full_bytes = num_dest_samples / 4;
  auto num_samples_left = num_dest_samples % 4;
  auto num_bytes = num_full_bytes + (num_samples_left > 0 ? 1 : 0);
  std::vector<int32_t> idx(num_dest_sample);
  std::copy(sample_idx, sample_idx + num_dest_samples, idx.begin());
  plhs[0] = mxCreateNumericMatrix(num_bytes, num_snps, mxUINT8_CLASS, mxREAL);
  auto *dest_geno = (uint8_t *)mxGetData(plhs[0]);
  snplib::Keep(src_geno, dest_geno, num_src_samples, num_dest_samples, num_snps,
               idx);
}
