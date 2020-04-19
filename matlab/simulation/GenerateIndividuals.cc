#include "../../src/simulations.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto *af = mxGetPr(prhs[0]);
  auto num_snps = static_cast<size_t>(mxGetM(prhs[0]));
  auto num_samples = static_cast<size_t>(mxGetScalar(prhs[1]));
  auto num_full_bytes = num_samples / 4;
  auto num_lefts = num_samples % 4;
  auto num_bytes = num_full_bytes + num_lefts != 0 ? 1 : 0;
  plhs[0] = mxCreateNumericMatrix(num_bytes, num_snps, mxUINT8_CLASS, mxREAL);
  auto *geno = reinterpret_cast<uint8_t *>(mxGetData(plhs[0]));
  snplib::GenerateIndividuals(af, num_samples, num_snps, geno);
}