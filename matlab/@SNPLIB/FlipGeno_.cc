#include "../../src//data_manage.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto *geno = (uint8_t *)mxGetData(prhs[0]);
  size_t num_samples = static_cast<size_t>(mxGetScalar(prhs[1]));
  size_t num_snps = static_cast<size_t>(mxGetM(prhs[2]));
  int32_t *snp_idx = (int32_t *)mxGetData(prhs[2]);
  std::vector<int32_t> idx(num_snps);
  std::copy(snp_idx, snp_idx + num_snps, idx.begin());
  snplib::FlipGeno(geno, num_samples, idx);
}
