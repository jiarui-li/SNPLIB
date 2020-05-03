#include "ugrm.h"

namespace {
const uint64_t kMask1 = 0x5555555555555555ull;  // 0x5555=b0101010101010101
const uint64_t kMask2 = 0xAAAAAAAAAAAAAAAAull;  // 0xAAAA=b1010101010101010
int64_t CalcUGRM(const __m128i *geno1, const __m128i *mask1,
                 const __m128i *geno2, const __m128i *mask2) noexcept {
  int64_t count_same = 0;
  int64_t count_diff = 0;
  for (size_t i = 0; i < 16; ++i) {
    __m128i g =
        _mm_xor_si128(_mm_load_si128(geno1 + i), _mm_load_si128(geno2 + i));
    __m128i m =
        _mm_and_si128(_mm_load_si128(mask1 + i), _mm_load_si128(mask2 + i));
    __m128i g_same = _mm_andnot_si128(g, m);
    __m128i g_diff = _mm_and_si128(g, m);
    count_same += _mm_popcnt_u64(_mm_cvtsi128_si64(g_same));
    count_diff += _mm_popcnt_u64(_mm_cvtsi128_si64(g_diff));
    g_same = _mm_srli_si128(g_same, 8);
    g_diff = _mm_srli_si128(g_diff, 8);
    count_same += _mm_popcnt_u64(_mm_cvtsi128_si64(g_same));
    count_diff += _mm_popcnt_u64(_mm_cvtsi128_si64(g_diff));
  }
  return count_same - count_diff;
}
void Mask(const uint64_t *geno64, size_t num_samples, uint64_t *mask64) {
  for (size_t i = 0; i < 32 * num_samples; ++i) {
    uint64_t tmp_geno = geno64[i] ^ kMask1;
    tmp_geno = (tmp_geno | (tmp_geno >> 1)) & kMask1;
    auto mask_1 = tmp_geno * 3ull;
    tmp_geno = geno64[i] ^ kMask2;
    tmp_geno = (tmp_geno | (tmp_geno << 1)) >> 1;
    auto mask_2 = (tmp_geno & kMask1) * 3ull;
    mask64[i] = mask_1 & mask_2;
  }
}
void UpdateMatrix(const uint64_t *geno64, const uint64_t *mask64,
                  size_t num_samples, int64_t *matrix) {
  for (size_t i = 0; i < num_samples; ++i) {
    auto *g1 = reinterpret_cast<const __m128i *>(geno64 + 32 * i);
    auto *m1 = reinterpret_cast<const __m128i *>(mask64 + 32 * i);
    for (size_t j = i; j < num_samples; ++j) {
      auto *g2 = reinterpret_cast<const __m128i *>(geno64 + 32 * j);
      auto *m2 = reinterpret_cast<const __m128i *>(mask64 + 32 * j);
      matrix[i * num_samples + j] += CalcUGRM(g1, m1, g2, m2);
    }
  }
}
}  // namespace

namespace snplib {
void CalcUGRMMatrixThread(const uint8_t *geno, size_t num_samples,
                          size_t num_snps, int64_t *matrix) {
  auto num_blocks = num_snps / 1024u;
  auto num_snps_left = num_snps % 1024u;
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + num_samples_left > 0 ? 1 : 0;
  auto *geno64 = new uint64_t[32 * num_samples];
  auto *mask64 = new uint64_t[32 * num_samples];
  std::fill(matrix, matrix + num_samples * num_samples, 0ull);
  for (size_t i = 0; i < num_blocks; ++i) {
    for (size_t j = 0; j < 32; ++j) {
      TransposeGeno(geno, num_samples, 32, j, geno64);
      geno += 32 * num_bytes;
    }
    Mask(geno64, num_samples, mask64);
    UpdateMatrix(geno64, mask64, num_samples, matrix);
  }
  if (num_snps_left > 0) {
    num_blocks = num_snps_left / 32;
    std::fill(geno64, geno64 + 32 * num_samples, kMask1);
    for (size_t j = 0; j < num_blocks; ++j) {
      TransposeGeno(geno, num_samples, 32, j, geno64);
      geno += 32 * num_bytes;
    }
    num_snps_left %= 32;
    if (num_snps_left > 0u) {
      TransposeGeno(geno, num_samples, num_snps_left, num_blocks, geno64);
    }
    Mask(geno64, num_samples, mask64);
    UpdateMatrix(geno64, mask64, num_samples, matrix);
  }
  delete[] geno64;
  delete[] mask64;
}
void CalcUGRMMatrix(const uint8_t *geno, size_t num_samples, size_t num_snps,
                    double *matrix, size_t num_threads) {
  std::vector<std::thread> workers;
  auto *matrices = new int64_t[num_samples * num_samples * num_threads];
  auto num_snps_job = num_snps / num_threads + 1;
  auto num_snps_left = num_snps % num_threads;
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + (num_samples_left > 0 ? 1 : 0);
  for (size_t i = 0; i < num_snps_left; ++i) {
    workers.emplace_back(CalcUGRMMatrixThread, geno, num_samples, num_snps_job,
                         matrices + i * num_samples * num_samples);
    geno += num_snps_job * num_bytes;
  }
  --num_snps_job;
  for (size_t i = num_snps_left; i < num_threads; ++i) {
    workers.emplace_back(CalcUGRMMatrixThread, geno, num_samples, num_snps_job,
                         matrices + i * num_samples * num_samples);
    geno += num_snps_job * num_bytes;
    geno += num_snps_job * num_bytes;
  }
  for (auto &&iter : workers) {
    iter.join();
  }
  auto *matrix_u = matrices;
  for (size_t k = 1; k < num_threads; ++k) {
    auto *tmp_m = matrices + k * num_samples * num_samples;
    for (size_t i = 0; i < num_samples; ++i) {
      for (size_t j = i; j < num_samples; ++j) {
        matrix_u[i * num_samples + j] += tmp_m[i * num_samples + j];
      }
    }
  }
  for (size_t i = 0; i < num_samples; ++i) {
    for (size_t j = i; j < num_samples; ++j) {
      matrix[i * num_samples + j] =
          static_cast<double>(matrix_u[i * num_samples + j]) / num_snps / 2.0;
      matrix[j * num_samples + i] = matrix[i * num_samples + j];
    }
  }
  delete[] matrices;
}
}  // namespace snplib