#ifndef SNPLIB_SRC_TRANSPOSE_GENO_H_
#define SNPLIB_SRC_TRANSPOSE_GENO_H_

#include <array>

namespace snplib {
void TransposeGeno(const uint8_t *geno, size_t num_samples, size_t num_snps,
                   size_t index, uint64_t *geno64);
void TransposeGeno(const uint8_t *geno, size_t num_samples, size_t num_snps,
                   uint64_t *geno64);
}  // namespace snplib

#endif  // SNPLIB_SRC_TRANSPOSE_GENO_H_
