#include "data_manage.h"

namespace {
uint8_t get(const uint8_t *geno, size_t index) {
  auto i = index / 4;
  auto s = index % 4;
  return (geno[i] >> (2 * s)) & 3u;
}
}  // namespace

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
              std::vector<int32_t> &index) {
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + num_samples_left > 0 ? 1 : 0;
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
          size_t num_dest_samples, size_t num_snps,
          std::vector<int32_t> &index) {
  auto num_full_bytes = num_src_samples / 4;
  auto num_samples_left = num_src_samples % 4;
  auto num_src_bytes = num_full_bytes + num_samples_left > 0 ? 1 : 0;
  num_full_bytes = num_dest_samples / 4;
  num_samples_left = num_dest_samples % 4;
  auto num_dest_bytes = num_full_bytes + num_samples_left > 0 ? 1 : 0;
  for (size_t i = 0; i < num_snps; ++i) {
    auto *tmp_dest_g = dest_geno + i * num_dest_bytes;
    auto *tmp_src_g = src_geno + i * num_src_bytes;
    for (size_t j = 0; j < num_full_bytes; ++j) {
      auto t = get(tmp_src_g, index[4 * j]);
      t += get(tmp_src_g, index[4 * j + 1]) << 2;
      t += get(tmp_src_g, index[4 * j + 2]) << 4;
      t += get(tmp_src_g, index[4 * j + 3]) << 6;
      tmp_dest_g[j] = t;
    }
    if (num_samples_left > 0) {
      uint8_t t = 0;
      for (size_t j = 0; j < num_samples_left; ++j) {
        t += get(tmp_src_g, index[4 * num_full_bytes + j]) << (2 * j);
      }
      tmp_dest_g[num_full_bytes] = t;
    }
  }
}
}  // namespace snplib