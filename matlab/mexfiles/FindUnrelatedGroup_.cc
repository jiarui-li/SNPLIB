#include "../../src/king.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto *matrix = mxGetPr(prhs[0]);
  auto num_samples = static_cast<size_t>(mxGetM(prhs[0]));
  auto threshold = mxGetScalar(prhs[1]);
  auto index = snplib::FindUnrelatedGroup(matrix, num_samples, threshold);
  auto length = index.size();
  plhs[0] = mxCreateDoubleMatrix(length, 1, mxREAL);
  auto *list = mxGetPr(plhs[0]);
  auto iter = index.begin();
  for (size_t i = 0; i < length; ++i) {
    list[i] = static_cast<double>(*iter + 1);
    iter++;
  }
}