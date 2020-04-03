#include "ibs.h"

namespace {
const uint64_t kMask = 0x5555555555555555ull;  // 0x5555=b0101010101010101
uint64_t CalcIBS(const __m128i *geno1, const __m128i *mask1,
                 const __m128i *geno2, const __m128i *mask2) {
  uint64_t count = 0;
  for (size_t i = 0; i < 16; ++i) {
    __m128i g =
        _mm_xor_si128(_mm_load_si128(geno1 + i), _mm_load_si128(geno2 + i));
    __m128i m =
        _mm_and_si128(_mm_load_si128(mask1 + i), _mm_load_si128(mask2 + i));
    g = _mm_andnot_si128(g, m);
    count += _mm_popcnt_u64(_mm_cvtsi128_si64(g));
    g = _mm_srli_si128(g, 8);
    count += _mm_popcnt_u64(_mm_cvtsi128_si64(g));
  }
  return count;
}
void Mask(const uint64_t *geno64, size_t num_samples, uint64_t *mask64) {
  for (size_t i = 0; i < 32 * num_samples; ++i) {
    uint64_t tmp_geno = geno64[i] ^ kMask;
    tmp_geno = (tmp_geno | (tmp_geno >> 1)) & kMask;
    mask64[i] = tmp_geno * 3ull;
  }
}
void UpdateMatrix(const uint64_t *geno64, const uint64_t *mask64,
                  size_t num_samples, uint64_t *matrix) {
  for (size_t i = 0; i < num_samples; ++i) {
    auto *g1 = reinterpret_cast<const __m128i *>(geno64 + 32 * i);
    auto *m1 = reinterpret_cast<const __m128i *>(mask64 + 32 * i);
    for (size_t j = i; j < num_samples; ++j) {
      auto *g2 = reinterpret_cast<const __m128i *>(geno64 + 32 * j);
      auto *m2 = reinterpret_cast<const __m128i *>(mask64 + 32 * j);
      matrix[i * num_samples + j] += CalcIBS(g1, m1, g2, m2);
    }
  }
}
void UpdateConnect(const uint64_t *src_geno64, const uint64_t *src_mask64,
                   const uint64_t *dest_geno64, const uint64_t *dest_mask64,
                   size_t num_src_samples, size_t num_dest_samples,
                   uint64_t *connect) {
  for (size_t i = 0; i < num_dest_samples; ++i) {
    auto *g1 = reinterpret_cast<const __m128i *>(dest_geno64 + 32 * i);
    auto *m1 = reinterpret_cast<const __m128i *>(dest_mask64 + 32 * i);
    uint64_t result = 0;
    for (size_t j = 0; j < num_src_samples; ++j) {
      auto *g2 = reinterpret_cast<const __m128i *>(src_geno64 + 32 * j);
      auto *m2 = reinterpret_cast<const __m128i *>(src_mask64 + 32 * j);
      result += CalcIBS(g1, m1, g2, m2);
    }
    connect[i] += result;
  }
}
void CalcIBSMatrixThread(uint8_t *geno, size_t num_samples, size_t num_snps,
                         uint64_t *matrix) {
  snplib::SNP snp(geno, num_samples);
  auto num_blocks = num_snps / 1024u;
  auto num_snps_left = num_snps % 1024u;
  auto *geno64 = new uint64_t[32 * num_samples];
  auto *mask64 = new uint64_t[32 * num_samples];
  std::fill(matrix, matrix + num_samples * num_samples, 0ull);
  for (size_t i = 0; i < num_blocks; ++i) {
    for (size_t j = 0; j < 32; ++j) {
      snplib::TransposeUGeno(snp, 32, j, geno64);
      snp += 32;
    }
    Mask(geno64, num_samples, mask64);
    UpdateMatrix(geno64, mask64, num_samples, matrix);
  }
  if (num_snps_left > 0) {
    num_blocks = num_snps_left / 32;
    std::fill(geno64, geno64 + 32 * num_samples, kMask);
    for (size_t j = 0; j < num_blocks; ++j) {
      snplib::TransposeUGeno(snp, 32, j, geno64);
      snp += 32;
    }
    num_snps_left %= 32;
    if (num_snps_left > 0u) {
      snplib::TransposeUGeno(snp, num_snps_left, num_blocks, geno64);
    }
    Mask(geno64, num_samples, mask64);
    UpdateMatrix(geno64, mask64, num_samples, matrix);
  }
  delete[] geno64;
  delete[] mask64;
}
void CalcIBSConnectThread(uint8_t *src_geno, size_t num_src_samples,
                          uint8_t *dest_geno, size_t num_dest_samples,
                          size_t num_snps, uint64_t *connection) {
  snplib::SNP src_snp(src_geno, num_src_samples);
  snplib::SNP dest_snp(dest_geno, num_dest_samples);
  auto num_blocks = num_snps / 1024u;
  auto num_snps_left = num_snps % 1024u;
  auto *src_geno64 = new uint64_t[32 * num_src_samples];
  auto *src_mask64 = new uint64_t[32 * num_src_samples];
  auto *dest_geno64 = new uint64_t[32 * num_dest_samples];
  auto *dest_mask64 = new uint64_t[32 * num_dest_samples];
  std::fill(connection, connection + num_dest_samples, 0ull);
  for (size_t i = 0; i < num_blocks; ++i) {
    for (size_t j = 0; j < 32; ++j) {
      snplib::TransposeUGeno(src_snp, 32, j, src_geno64);
      snplib::TransposeUGeno(dest_snp, 32, j, dest_geno64);
      src_snp += 32;
      dest_snp += 32;
    }
    Mask(src_geno64, num_src_samples, src_mask64);
    Mask(dest_geno64, num_dest_samples, dest_mask64);
    UpdateConnect(src_geno64, src_mask64, dest_geno64, dest_mask64,
                  num_src_samples, num_dest_samples, connection);
  }
  if (num_snps_left > 0u) {
    std::fill(src_geno64, src_geno64 + 32 * num_src_samples, kMask);
    std::fill(dest_geno64, dest_geno64 + 32 * num_dest_samples, kMask);
    num_blocks = num_snps_left / 32;
    for (size_t j = 0; j < num_blocks; ++j) {
      snplib::TransposeUGeno(src_snp, 32, j, src_geno64);
      snplib::TransposeUGeno(dest_snp, 32, j, dest_geno64);
      src_snp += 32;
      dest_snp += 32;
    }
    num_snps_left %= 32;
    if (num_snps_left > 0u) {
      snplib::TransposeUGeno(src_snp, num_snps_left, num_blocks, src_geno64);
      snplib::TransposeUGeno(dest_snp, num_snps_left, num_blocks, dest_geno64);
    }
    Mask(src_geno64, num_src_samples, src_mask64);
    Mask(dest_geno64, num_dest_samples, dest_mask64);
    UpdateConnect(src_geno64, src_mask64, dest_geno64, dest_mask64,
                  num_src_samples, num_dest_samples, connection);
  }
  delete[] src_geno64;
  delete[] src_mask64;
  delete[] dest_geno64;
  delete[] dest_mask64;
}
}  // namespace
void CalcIBSMatrix(uint8_t *geno, size_t num_samples, size_t num_snps,
                   double *matrix, size_t num_threads) {
  std::vector<std::thread> workers(num_threads);
  auto *matrices = new uint64_t[num_samples * num_samples * num_threads];
  auto num_snps_job = num_snps / num_threads + 1;
  auto num_snps_left = num_snps % num_threads;
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + (num_samples_left > 0 ? 1 : 0);
  for (size_t i = 0; i < num_snps_left; ++i) {
    workers[i] =
        std::thread(CalcIBSMatrixThread, geno, num_samples, num_snps_job,
                    matrices + i * num_samples * num_samples);
    geno += num_snps_job * num_bytes;
  }
  --num_snps_job;
  for (size_t i = num_snps_left; i < num_threads; ++i) {
    workers[i] =
        std::thread(CalcIBSMatrixThread, geno, num_samples, num_snps_job,
                    matrices + i * num_samples * num_samples);
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
void CalcIBSConnection(uint8_t *src_geno, size_t num_src_samples,
                       uint8_t *dest_geno, size_t num_dest_samples,
                       size_t num_snps, double *connection,
                       size_t num_threads) {
  std::vector<std::thread> workers(num_threads);
  auto *connects = new uint64_t[num_dest_samples * num_threads];
  auto num_src_full_bytes = num_src_samples / 4;
  auto num_src_samples_left = num_src_samples % 4;
  auto num_src_bytes = num_src_full_bytes + (num_src_samples_left > 0 ? 1 : 0);
  auto num_dest_full_bytes = num_dest_samples / 4;
  auto num_dest_samples_left = num_dest_samples % 4;
  auto num_dest_bytes =
      num_dest_full_bytes + (num_dest_samples_left > 0 ? 1 : 0);
  auto num_snps_job = num_snps / num_threads + 1;
  auto num_snps_left = num_snps % num_threads;
  for (size_t i = 0; i < num_snps_left; ++i) {
    workers[i] = std::thread(CalcIBSConnectThread, src_geno, num_src_samples,
                             dest_geno, num_dest_samples, num_snps_job,
                             connects + i * num_dest_samples);
    src_geno += num_snps_job * num_src_bytes;
    dest_geno += num_snps_job * num_dest_bytes;
  }
  --num_snps_job;
  for (size_t i = num_snps_left; i < num_threads; ++i) {
    workers[i] = std::thread(CalcIBSConnectThread, src_geno, num_src_samples,
                             dest_geno, num_dest_samples, num_snps_job,
                             connects + i * num_dest_samples);
    src_geno += num_snps_job * num_src_bytes;
    dest_geno += num_snps_job * num_dest_bytes;
  }
  for (auto &&iter : workers) {
    iter.join();
  }
  auto *connect_u = connects;
  for (size_t i = 1; i < num_threads; ++i) {
    auto *tmp_c = connects + i * num_dest_samples;
    for (size_t j = 0; j < num_dest_samples; ++j) {
      connect_u[j] += tmp_c[j];
    }
  }
  for (size_t i = 0; i < num_dest_samples; ++i) {
    connection[i] = static_cast<double>(connect_u[i]) / num_snps / 2.0;
  }
  delete[] connects;
}