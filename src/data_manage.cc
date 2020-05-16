#include "data_manage.h"

namespace snplib {
void UnpackGeno(const uint8_t *geno, size_t num_samples, size_t num_snps,
                double *geno_d) {
  const std::array<double, 4> geno_table{0.0, 0.0, 1.0, 2.0};
  SNP snp(geno, num_samples);
  for (size_t i = 0; i < num_snps; ++i) {
    auto *snp_geno_d = geno_d + i * num_samples;
    snp.UnpackGeno(geno_table, snp_geno_d);
    snp += 1;
  }
}
void UnpackGRMGeno(const uint8_t *geno, const double *af, size_t num_samples,
                   size_t num_snps, double *geno_d) {
  std::array<double, 4> geno_table{0.0, 0.0, 0.0, 0.0};
  SNP snp(geno, num_samples);
  for (size_t i = 0; i < num_snps; ++i) {
    auto mu = 2.0 * af[i];
    auto std = std::sqrt(2.0 * af[i] * (1.0 - af[i]));
    geno_table[0] = -mu / std;
    geno_table[2] = (1.0 - mu) / std;
    geno_table[3] = (2.0 - mu) / std;
    auto *snp_geno_d = geno_d + i * num_samples;
    snp.UnpackGeno(geno_table, snp_geno_d);
    snp += 1;
  }
}
void UnpackUGeno(const uint8_t *geno, size_t num_samples, size_t num_snps,
                 double *geno_d) {
  const std::array<double, 4> geno_table{-1.0, 0.0, 0.0, 1.0};
  SNP snp(geno, num_samples);
  for (size_t i = 0; i < num_snps; ++i) {
    auto *snp_geno_d = geno_d + i * num_samples;
    snp.UnpackGeno(geno_table, snp_geno_d);
    snp += 1;
  }
}
void FlipGeno(uint8_t *geno, size_t num_samples, size_t num_snps,
              const int32_t *index) {
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + (num_samples_left > 0 ? 1 : 0);
  for (size_t i = 0; i < num_snps; ++i) {
    auto *tmp_g = geno + index[i] * num_bytes;
    for (size_t j = 0; j < num_bytes; ++j) {
      uint8_t g = ~tmp_g[j];
      uint8_t g1 = (g >> 1) & (uint8_t)0x55;
      uint8_t g2 = (g << 1) & (uint8_t)0xAA;
      tmp_g[j] = g1 + g2;
    }
  }
}
void Keep(const uint8_t *src_geno, uint8_t *dest_geno, size_t num_src_samples,
          size_t num_dest_samples, size_t num_snps, const int32_t *index) {
  SNP src_snp(src_geno, num_src_samples);
  auto num_full_bytes = num_dest_samples / 4;
  auto num_samples_left = num_dest_samples % 4;
  auto num_bytes = num_full_bytes + (num_samples_left > 0 ? 1 : 0);
  for (size_t l = 0; l < num_snps; ++l) {
    auto *tmp_g = dest_geno + l * num_bytes;
    for (size_t i = 0; i < num_full_bytes; ++i) {
      uint8_t t = src_snp(index[4 * i]);
      t += src_snp(index[4 * i + 1]) << 2;
      t += src_snp(index[4 * i + 2]) << 4;
      t += src_snp(index[4 * i + 3]) << 6;
      tmp_g[i] = t;
    }
    if (num_samples_left > 0) {
      uint8_t t = 0;
      for (size_t i = 0; i < num_samples_left; ++i) {
        t += src_snp(index[4 * num_full_bytes + i]) << (2 * i);
      }
      tmp_g[num_full_bytes] = t;
    }
    src_snp += 1;
  }
}
}  // namespace snplib