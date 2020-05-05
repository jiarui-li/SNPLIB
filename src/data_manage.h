#ifndef SNPLIB_SRC_DATA_MANAGE_H_
#define SNPLIB_SRC_DATA_MANAGE_H_

#include <cmath>
#include <vector>

#include "unpack_geno.h"

namespace snplib {
void UnpackGRMGeno(const uint8_t *geno, const double *af, size_t num_samples,
                   size_t num_snps, double *geno_d);
void UnpackUGeno(const uint8_t *geno, size_t num_samples, size_t num_snps,
                 double *geno_d);
void FlipGeno(uint8_t *geno, size_t num_samples, size_t num_snps,
              std::vector<int32_t> &index);
void Keep(const uint8_t *src_geno, uint8_t *dest_geno, size_t num_src_samples,
          size_t num_dest_samples, size_t num_snps,
          std::vector<int32_t> &index);
}  // namespace snplib

#endif  // SNPLIB_SRC_DATA_MANAGE_H_
