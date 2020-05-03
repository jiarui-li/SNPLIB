#ifndef SNPLIB_SRC_UNPACK_GENO_H_
#define SNPLIB_SRC_UNPACK_GENO_H_

#include <array>

namespace snplib {
template <class T>
void UnpackGeno(const uint8_t *geno, size_t num_samples,
                const std::array<T, 4> &geno_table,
                const std::array<T, 4> &mask_table, T *geno_d, T *mask_d) {
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + num_samples_left > 0 ? 1 : 0;
  for (size_t i = 0; i < num_full_bytes; ++i) {
    auto t = geno[i];
    geno_d[4 * i] = geno_table[t & 3u];
    mask_d[4 * i] = mask_table[t & 3u];
    t >>= 2;
    geno_d[4 * i + 1] = geno_table[t & 3u];
    mask_d[4 * i + 1] = mask_table[t & 3u];
    t >>= 2;
    geno_d[4 * i + 2] = geno_table[t & 3u];
    mask_d[4 * i + 2] = mask_table[t & 3u];
    t >>= 2;
    geno_d[4 * i + 3] = geno_table[t & 3u];
    mask_d[4 * i + 3] = mask_table[t & 3u];
  }
  if (num_samples_left > 0) {
    auto t = geno[num_full_bytes];
    for (size_t i = 0; i < num_samples_left; ++i) {
      geno_d[4 * num_full_bytes + 2] = geno_table[t & 3u];
      mask_d[4 * num_full_bytes + 2] = mask_table[t & 3u];
      t >>= 2;
    }
  }
}

template <class T>
void UnpackGeno(const uint8_t *geno, size_t num_samples,
                const std::array<T, 4> &geno_table, T *geno_d) {
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + num_samples_left > 0 ? 1 : 0;
  for (size_t i = 0; i < num_full_bytes; ++i) {
    auto t = geno[i];
    geno_d[4 * i] = geno_table[t & 3u];
    t >>= 2;
    geno_d[4 * i + 1] = geno_table[t & 3u];
    t >>= 2;
    geno_d[4 * i + 2] = geno_table[t & 3u];
    t >>= 2;
    geno_d[4 * i + 3] = geno_table[t & 3u];
  }
  if (num_samples_left > 0) {
    auto t = geno[num_full_bytes];
    for (size_t i = 0; i < num_samples_left; ++i) {
      geno_d[4 * num_full_bytes + 2] = geno_table[t & 3u];
      t >>= 2;
    }
  }
}
}  // namespace snplib

#endif  // SNPLIB_SRC_UNPACK_GENO_H_
