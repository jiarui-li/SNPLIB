#include "snp.h"

namespace {
const double kGenoTable[4] = {0.0, 0.0, 1.0, 2.0};
const double kMaskTable[4] = {1.0, 0.0, 1.0, 1.0};
const uint64_t kTransUG = 0x0303030303030303ull;  // 0303 = 0000001100000011
const uint64_t kTransSG = 0x0003000300030003ull;  // 0003 = 0000000000000011
uint64_t ConvertUG64(const std::array<uint64_t, 4> &g64) {
  uint64_t geno64;
  geno64 = g64[3] & kTransUG;
  geno64 <<= 2;
  geno64 |= g64[2] & kTransUG;
  geno64 <<= 2;
  geno64 |= g64[1] & kTransUG;
  geno64 <<= 2;
  geno64 |= g64[0] & kTransUG;
  return geno64;
}
uint64_t ConvertSG64(const std::array<uint64_t, 5> &g64) {
  uint64_t geno64;
  geno64 = g64[4] & kTransSG;
  geno64 <<= 3;
  geno64 |= g64[3] & kTransSG;
  geno64 <<= 3;
  geno64 |= g64[2] & kTransSG;
  geno64 <<= 3;
  geno64 |= g64[1] & kTransSG;
  geno64 <<= 3;
  geno64 |= g64[0] & kTransSG;
  return geno64;
}
}  // namespace

namespace snplib {
void SNP::ConvertUGeno(size_t num_snps, size_t idx,
                       std::array<uint64_t, 4> &g64) {
  uint64_t g8[32] = {
      0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull,
      0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull,
      0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull,
      0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull};
  for (size_t i = 0; i < num_snps; ++i) {
    auto *tmp_g = geno_ + i * num_bytes_;
    g8[i] = static_cast<uint64_t>(tmp_g[idx]);
  }
  for (size_t i = 0; i < 4; ++i) {
    g64[i] = g8[i] + (g8[4 + i] << 8) + (g8[8 + i] << 16) + (g8[12 + i] << 24) +
             (g8[16 + i] << 32) + (g8[20 + i] << 40) + (g8[24 + i] << 48) +
             (g8[28 + i] << 56);
  }
}
void SNP::ConvertSGeno(size_t num_snps, size_t idx,
                       std::array<uint64_t, 5> &g64) {
  uint64_t g8[20] = {0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull,
                     0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull,
                     0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull,
                     0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull};
  for (size_t i = 0; i < num_snps; ++i) {
    auto *tmp_g = geno_ + i * num_bytes_;
    g8[i] = static_cast<uint64_t>(tmp_g[idx]);
  }
  for (size_t i = 0; i < 5; ++i) {
    g64[i] =
        g8[i] + (g8[5 + i] << 16) + (g8[10 + i] << 32) + (g8[15 + i] << 48);
  }
}
void SNP::Flip() {
  for (size_t i = 0; i < num_bytes_; ++i) {
    uint8_t g = ~geno_[i];
    uint8_t g1 = (g >> 1) & (uint8_t)0x55;
    uint8_t g2 = (g << 1) & (uint8_t)0xAA;
    geno_[i] = g1 + g2;
  }
}
void SNP::UnpackGeno(double *geno_d, double *mask_d) {
  for (size_t i = 0; i < num_full_bytes_; ++i) {
    auto t = geno_[i];
    geno_d[4 * i] = kGenoTable[t & 3u];
    mask_d[4 * i] = kMaskTable[t & 3u];
    t >>= 2;
    geno_d[4 * i + 1] = kGenoTable[t & 3u];
    mask_d[4 * i + 1] = kMaskTable[t & 3u];
    t >>= 2;
    geno_d[4 * i + 2] = kGenoTable[t & 3u];
    mask_d[4 * i + 2] = kMaskTable[t & 3u];
    t >>= 2;
    geno_d[4 * i + 3] = kGenoTable[t & 3u];
    mask_d[4 * i + 3] = kMaskTable[t & 3u];
  }
  if (num_samples_left_ > 0u) {
    auto t = geno_[num_full_bytes_];
    for (size_t i = 0; i < num_samples_left_; ++i) {
      geno_d[4 * num_full_bytes_ + i] = kGenoTable[t & 3u];
      mask_d[4 * num_full_bytes_ + i] = kMaskTable[t & 3u];
      t >>= 2;
    }
  }
}
void SNP::UnpackGeno(double *geno_d, const std::array<double, 4> &geno_table) {
  for (size_t i = 0; i < num_full_bytes_; ++i) {
    auto t = geno_[i];
    geno_d[4 * i] = geno_table[t & 3u];
    t >>= 2;
    geno_d[4 * i + 1] = geno_table[t & 3u];
    t >>= 2;
    geno_d[4 * i + 2] = geno_table[t & 3u];
    t >>= 2;
    geno_d[4 * i + 3] = geno_table[t & 3u];
  }
  if (num_samples_left_ > 0u) {
    auto t = geno_[num_full_bytes_];
    for (size_t i = 0; i < num_samples_left_; ++i) {
      geno_d[4 * num_full_bytes_ + i] = geno_table[t & 3u];
      t >>= 2;
    }
  }
}
void SNP::TransposeUGeno(size_t num_snps, size_t idx, uint64_t *geno64) {
  std::array<uint64_t, 4> g64;
  for (size_t i = 0; i < num_full_bytes_; ++i) {
    ConvertUGeno(num_snps, i, g64);
    geno64[32 * (4 * i) + idx] = ConvertUG64(g64);
    for (size_t j = 1; j < 4; ++j) {
      for (auto &iter : g64) {
        iter >>= 2;
      }
      geno64[32 * (4 * i + j) + idx] = ConvertUG64(g64);
    }
  }
  if (num_samples_left_ > 0) {
    ConvertUGeno(num_snps, num_full_bytes_, g64);
    geno64[32 * (4 * num_full_bytes_) + idx] = ConvertUG64(g64);
    for (size_t j = 1; j < num_samples_left_; ++j) {
      for (auto &iter : g64) {
        iter >>= 2;
      }
      geno64[32 * (4 * num_full_bytes_ + j) + idx] = ConvertUG64(g64);
    }
  }
}
void SNP::TransposeSGeno(size_t num_snps, uint64_t *geno64) {
  std::array<uint64_t, 5> g64;
  for (size_t i = 0; i < num_full_bytes_; ++i) {
    ConvertSGeno(num_snps, i, g64);
    geno64[4 * i] = ConvertSG64(g64);
    for (size_t j = 1; j < 4; ++j) {
      for (auto &iter : g64) {
        iter >>= 2;
      }
      geno64[4 * i + j] = ConvertSG64(g64);
    }
  }
  if (num_samples_left_ > 0) {
    ConvertSGeno(num_snps, num_full_bytes_, g64);
    geno64[4 * num_full_bytes_] = ConvertSG64(g64);
    for (size_t j = 1; j < num_samples_left_; ++j) {
      for (auto &iter : g64) {
        iter >>= 2;
      }
      geno64[4 * num_full_bytes_ + j] = ConvertSG64(g64);
    }
  }
}

}  // namespace snplib