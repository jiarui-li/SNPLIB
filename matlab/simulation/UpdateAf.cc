#include "../../src/simulations.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto *aaf = mxGetPr(prhs[0]);
  auto num_snps = static_cast<size_t>(mxGetM(prhs[0]));
  auto num_pops = static_cast<size_t>(mxGetScalar(prhs[1]));
  auto num_generations = static_cast<size_t>(mxGetScalar(prhs[2]));
  auto Ne = static_cast<size_t>(mxGetScalar(prhs[3]));
  plhs[0] = mxCreateDoubleMatrix(num_snps, num_pops, mxREAL);
  auto *af = mxGetPr(plhs[0]);
  snplib::UpdateAf(aaf, num_pops, num_snps, num_generations, Ne, af);
}
