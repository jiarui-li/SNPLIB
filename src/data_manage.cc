#include "data_manage.h"

namespace snplib {
void FlipGeno(uint8_t *geno, size_t num_samples,
              const std::vector<int32_t> &idx) {
  SNP snp(geno, num_samples);
  for (auto iter : idx) {
    auto tmp_s = snp;
    tmp_s += iter;
    tmp_s.Flip();
  }
}
void Keep(uint8_t *src_geno, uint8_t *dest_geno, size_t num_src_samples,
          size_t num_dest_samples, size_t num_snps,
          const std::vector<int32_t> &idx) {
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
void UnpackGRMGeno(uint8_t *geno, const double *af, size_t num_samples,
                   size_t num_snps, double *grm_geno) {
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + (num_samples_left > 0 ? 1 : 0);
  for (size_t i = 0; i < num_snps; ++i) {
    auto *snp_geno = geno + i * num_bytes;
    auto *snp_geno_d = grm_geno + i * num_samples;
    auto standard_deviation = std::sqrt(2.0 * af[i] * (1.0 - af[i]));
    auto mu = 2.0 * af[i];
    double geno_table[4] = {-mu / standard_deviation, 0.0,
                            (1.0 - mu) / standard_deviation,
                            (2.0 - mu) / standard_deviation};
    for (size_t j = 0; j < num_full_bytes; ++j) {
      auto t = snp_geno[j];
      snp_geno_d[4 * j] = geno_table[t & 3u];
      t >>= 2;
      snp_geno_d[4 * j + 1] = geno_table[t & 3u];
      t >>= 2;
      snp_geno_d[4 * j + 2] = geno_table[t & 3u];
      t >>= 2;
      snp_geno_d[4 * j + 3] = geno_table[t & 3u];
    }
    if (num_samples_left > 0u) {
      auto t = snp_geno[num_full_bytes];
      for (size_t j = 0; j < num_samples_left; ++j) {
        snp_geno_d[4 * num_full_bytes + j] = geno_table[t & 3u];
        t >>= 2;
      }
    }
  }
}
void UnpackUGeno(uint8_t *geno, size_t num_samples, size_t num_snps,
                 double *geno_d) {
  const double geno_table[4] = {-1.0, 0.0, 0.0, 1.0};
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + (num_samples_left > 0 ? 1 : 0);
  for (size_t i = 0; i < num_snps; ++i) {
    auto *snp_geno = geno + i * num_bytes;
    auto *snp_geno_d = geno_d + i * num_samples;
    for (size_t j = 0; j < num_full_bytes; ++j) {
      auto t = snp_geno[j];
      snp_geno_d[4 * j] = geno_table[t & 3u];
      t >>= 2;
      snp_geno_d[4 * j + 1] = geno_table[t & 3u];
      t >>= 2;
      snp_geno_d[4 * j + 2] = geno_table[t & 3u];
      t >>= 2;
      snp_geno_d[4 * j + 3] = geno_table[t & 3u];
    }
    if (num_samples_left > 0u) {
      auto t = snp_geno[num_full_bytes];
      for (size_t j = 0; j < num_samples_left; ++j) {
        snp_geno_d[4 * num_full_bytes + j] = geno_table[t & 3u];
        t >>= 2;
      }
    }
  }
}
}  // namespace snplib