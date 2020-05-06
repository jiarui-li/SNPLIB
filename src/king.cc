#include "king.h"

namespace {
const uint64_t kMask1 = 0x5555555555555555ull;  // 0x5555=b0101010101010101
const uint64_t kMask2 = 0xAAAAAAAAAAAAAAAAull;  // 0xAAAA=b1010101010101010
int64_t CalcKING(const uint64_t *geno1, const uint64_t *geno2) {
  int64_t num_hetero = 0;
  int64_t num_homo = 0;
  for (size_t i = 0; i < 32; ++i) {
    auto mask_1 = geno1[i] ^ kMask1;
    mask_1 = (mask_1 | (mask_1 >> 1)) & kMask1;
    mask_1 *= 3ull;
    auto mask_2 = geno2[i] ^ kMask1;
    mask_2 = (mask_2 | (mask_2 >> 1)) & kMask1;
    mask_2 *= 3ull;
    auto hetero_1 = geno1[i] ^ kMask2;
    hetero_1 = (hetero_1 | (hetero_1 << 1)) >> 1;
    hetero_1 = (hetero_1 & kMask1) * 3ull;
    auto hetero_2 = geno2[i] ^ kMask2;
    hetero_2 = (hetero_2 | (hetero_2 << 1)) >> 1;
    hetero_2 = (hetero_2 & kMask1) * 3ull;
    num_hetero += _mm_popcnt_u64(~(hetero_1 | hetero_2));
    auto geno = geno1[i] ^ geno2[i];
    geno &= mask_1 & mask_2 & hetero_1 & hetero_2;
    num_homo += _mm_popcnt_u64(geno);
  }
  return num_hetero - 2 * num_homo;
}
uint64_t CalcHetero(const uint64_t *geno) {
  uint64_t num_hetero = 0;
  for (size_t i = 0; i < 32; ++i) {
    auto mask = geno[i] ^ kMask1;
    mask = (mask | (mask >> 1)) & kMask1;
    mask *= 3ull;
    auto hetero = geno[i] ^ kMask2;
    hetero = (hetero | (hetero << 1)) >> 1;
    hetero = (hetero & mask) * 3ull;
    num_hetero += _mm_popcnt_u64(~hetero);
  }
  return num_hetero;
}
void UpdateMatrix(const uint64_t *geno64, size_t num_samples, int64_t *matrix) {
  for (size_t i = 0; i < num_samples; ++i) {
    auto *g1 = geno64 + 32 * i;
    for (size_t j = i + 1; j < num_samples; ++j) {
      auto *g2 = geno64 + 32 * j;
      matrix[i * num_samples + j] += CalcKING(g1, g2);
    }
  }
}
void UpdateHetero(const uint64_t *geno64, size_t num_samples,
                  uint64_t *vector) {
  for (size_t i = 0; i < num_samples; ++i) {
    auto *g1 = geno64 + 32 * i;
    vector[i] += CalcHetero(g1);
  }
}
void CalcKINGThread(const uint8_t *geno, size_t num_samples, size_t num_snps,
                    int64_t *matrix, uint64_t *vector) {
  snplib::SNP snp(geno, num_samples);
  auto num_blocks = num_snps / 1024u;
  auto num_snps_left = num_snps % 1024u;
  auto *geno64 = new uint64_t[32 * num_samples];
  std::fill(matrix, matrix + num_samples * num_samples, 0ull);
  std::fill(vector, vector + num_samples, 0ull);
  for (size_t i = 0; i < num_blocks; ++i) {
    for (size_t j = 0; j < 32; ++j) {
      snp.TransposeGeno(32, j, geno64);
      snp += 32;
    }
    UpdateHetero(geno64, num_samples, vector);
    UpdateMatrix(geno64, num_samples, matrix);
  }
  if (num_snps_left > 0) {
    num_blocks = num_snps_left / 32;
    std::fill(geno64, geno64 + 32 * num_samples, kMask1);
    for (size_t j = 0; j < num_blocks; ++j) {
      snp.TransposeGeno(32, j, geno64);
      snp += 32;
    }
    num_snps_left %= 32;
    if (num_snps_left > 0u) {
      snp.TransposeGeno(num_snps_left, num_blocks, geno64);
    }
    UpdateHetero(geno64, num_samples, vector);
    UpdateMatrix(geno64, num_samples, matrix);
  }
  delete[] geno64;
}
}  // namespace

namespace snplib {
void CalcKINGMatrix(const uint8_t *geno, size_t num_samples, size_t num_snps,
                    double *matrix, size_t num_threads) {
  std::vector<std::thread> workers;
  auto *matrices = new int64_t[num_samples * num_samples * num_threads];
  auto *vectors = new uint64_t[num_samples * num_threads];
  auto num_snps_job = num_snps / num_threads + 1;
  auto num_snps_left = num_snps % num_threads;
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + (num_samples_left > 0 ? 1 : 0);
  for (size_t i = 0; i < num_snps_left; ++i) {
    workers.emplace_back(CalcKINGThread, geno, num_samples, num_snps_job,
                         matrices + i * num_samples * num_samples,
                         vectors + i * num_samples);
    geno += num_snps_job * num_bytes;
  }
  --num_snps_job;
  for (size_t i = num_snps_left; i < num_threads; ++i) {
    workers.emplace_back(CalcKINGThread, geno, num_samples, num_snps_job,
                         matrices + i * num_samples * num_samples,
                         vectors + i * num_samples);
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
  auto *vector_u = vectors;
  for (size_t k = 1; k < num_threads; ++k) {
    auto *tmp_v = vectors + k * num_samples;
    for (size_t i = 0; i < num_samples; ++i) {
      vector_u[i] += tmp_v[i];
    }
  }
  for (size_t i = 0; i < num_samples; ++i) {
    for (size_t j = i + 1; j < num_samples; ++j) {
      matrix[i * num_samples + j] =
          static_cast<double>(matrix_u[i * num_samples + j]) /
          (vector_u[i] + vector_u[j]);
      matrix[j * num_samples + i] = matrix[i * num_samples + j];
    }
  }
  delete[] matrices;
  delete[] vectors;
}
}  // namespace snplib