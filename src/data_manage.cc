#include "data_manage.h"

namespace snplib {
void FlipGeno(uint8_t *geno, size_t num_samples, size_t num_snps,
              const int32_t *idx) {
  SNP snp(geno, num_samples);
  for (size_t i = 0; i < num_snps; ++i) {
    auto tmp_s = snp;
    tmp_s += idx[i];
    tmp_s.Flip();
  }
}
void Keep(uint8_t *src_geno, uint8_t *dest_geno, size_t num_src_samples,
          size_t num_dest_samples, size_t num_snps, const int32_t *idx) {
  SNP src_snp(src_geno, num_src_samples);
  auto num_full_bytes = num_dest_samples / 4;
  auto num_samples_left = num_dest_samples % 4;
  auto num_bytes = num_full_bytes + (num_samples_left > 0 ? 1 : 0);
  for (size_t l = 0; l < num_snps; ++l) {
    auto *tmp_g = dest_geno + l * num_bytes;
    for (size_t i = 0; i < num_full_bytes; ++i) {
      uint8_t t = src_snp[idx[4 * i]];
      t += src_snp[idx[4 * i + 1]] << 2;
      t += src_snp[idx[4 * i + 2]] << 4;
      t += src_snp[idx[4 * i + 3]] << 6;
      tmp_g[i] = t;
    }
    if (num_samples_left > 0) {
      uint8_t t = 0;
      for (size_t i = 0; i < num_samples_left; ++i) {
        t += src_snp[idx[4 * num_full_bytes + i]] << (2 * i);
      }
      tmp_g[num_full_bytes] = t;
    }
    ++src_snp;
  }
}
}  // namespace snplib