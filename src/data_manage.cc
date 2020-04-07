#include "data_manage.h"

namespace snplib {

void FlipGeno(uint8_t *geno, size_t num_samples, size_t num_snps,
              const int32_t *idx) {}
void Keep(uint8_t *src_geno, uint8_t *dest_geno, size_t num_src_samples,
          size_t num_dest_samples, size_t num_snps, const int32_t *idx) {}
void UnpackGRMGeno(uint8_t *geno, const double *af, size_t num_samples,
                   size_t num_snps, double *grm_geno) {}
void UnpackUGeno(uint8_t *geno, size_t num_samples, size_t num_snps,
                 double *geno_d) {}
}  // namespace snplib