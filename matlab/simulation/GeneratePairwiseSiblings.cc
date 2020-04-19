#include "../../src/simulations.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto *parent_geno = reinterpret_cast<uint8_t *>(mxGetData(prhs[0]));
  auto num_snps = static_cast<size_t>(mxGetM(prhs[0]));
  auto num_families = static_cast<size_t>(mxGetScalar(prhs[1]));
  auto num_samples = 2 * num_families;
  auto num_full_bytes = num_samples / 2;
  auto num_lefts = num_samples % 4;
  auto num_bytes = num_full_bytes + num_lefts != 0 ? 1 : 0;
  plhs[0] = mxCreateNumericMatrix(num_bytes, num_snps, mxUINT8_CLASS, mxREAL);
  auto *geno = reinterpret_cast<uint8_t *>(mxGetData(plhs[0]));
  snplib::GeneratePairwiseSiblings(parent_geno, num_samples, num_snps, geno);
}