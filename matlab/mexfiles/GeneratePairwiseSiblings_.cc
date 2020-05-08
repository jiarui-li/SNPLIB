#include "../../src/simulations.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto *parent_geno = reinterpret_cast<uint8_t *>(mxGetData(prhs[0]));
  auto num_bytes = static_cast<size_t>(mxGetM(prhs[0]));
  auto num_snps = static_cast<size_t>(mxGetN(prhs[0]));
  auto num_samples = static_cast<size_t>(mxGetScalar(prhs[1]));
  auto num_families = num_samples / 2;
  plhs[0] = mxCreateNumericMatrix(num_bytes, num_snps, mxUINT8_CLASS, mxREAL);
  auto *siblings_geno = reinterpret_cast<uint8_t *>(mxGetData(plhs[0]));
  snplib::GeneratePairwiseSiblings(parent_geno, num_families, num_snps,
                                   siblings_geno);
}