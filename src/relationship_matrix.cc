#include "relationship_matrix.h"

namespace {
const uint64_t kMask1 = 0x5555555555555555ull;  // 0x5555=b0101010101010101
const uint64_t kMask2 = 0xAAAAAAAAAAAAAAAAull;  // 0xAAAA=b1010101010101010
const uint64_t kMask3 = 0x1249124912491249ull;  // 0x1249=b0001001001001001
}  // namespace

namespace ibs {
uint64_t CalcDistance(const __m128i *geno1, const __m128i *mask1,
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
    uint64_t tmp_geno = geno64[i] ^ kMask1;
    tmp_geno = (tmp_geno | (tmp_geno >> 1)) & kMask1;
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
      matrix[i * num_samples + j] += CalcDistance(g1, m1, g2, m2);
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
      result += CalcDistance(g1, m1, g2, m2);
    }
    connect[i] += result;
  }
}
}  // namespace ibs

namespace ugrm {
int64_t CalcDistance(const __m128i *geno1, const __m128i *mask1,
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
    tmp_geno = (tmp_geno | tmp_geno >> 1) & kMask1;
    auto mask_1 = tmp_geno * 3ull;
    tmp_geno = geno64[i] ^ kMask2;
    tmp_geno = (tmp_geno | tmp_geno << 1) >> 1;
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
      matrix[i * num_samples + j] += CalcDistance(g1, m1, g2, m2);
    }
  }
}
}  // namespace ugrm

namespace king {
int64_t CalcDistance(const uint64_t *geno1, const uint64_t *geno2) {
  int64_t num_hetero = 0;
  int64_t num_homo = 0;
  for (size_t i = 0; i < 32; ++i) {
    auto mask_1 = geno1[i] ^ kMask1;
    mask_1 = (mask_1 | mask_1 >> 1) & kMask1;
    mask_1 *= 3ull;
    auto mask_2 = geno2[i] ^ kMask1;
    mask_2 = (mask_2 | mask_2 >> 1) & kMask1;
    mask_2 *= 3ull;
    auto hetero_1 = geno1[i] ^ kMask2;
    hetero_1 = (hetero_1 | hetero_1 << 1) >> 1;
    hetero_1 = (hetero_1 & kMask1) * 3ull;
    auto hetero_2 = geno2[i] ^ kMask2;
    hetero_2 = (hetero_2 | hetero_2 << 1) >> 1;
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
    hetero = (hetero | hetero << 1) >> 1;
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
      matrix[i * num_samples + j] += CalcDistance(g1, g2);
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
}  // namespace king

namespace grm {
void CalcSnpCorr(double af, std::array<double, 8> &table) {
  table[0] = 2.0 * af / (1.0 - af);        // 00 00
  table[1] = 0.0;                          // 00 01
  table[2] = 1.0 / (1.0 - af) - 2;         // 00 10
  table[3] = -2.0;                         // 00 11
  table[4] = 0.5 / af / (1.0 - af) - 2.0;  // 10 10
  table[5] = 1.0 / af - 2.0;               // 10 11
  table[6] = 2.0 / af - 2.0;               // 11 11
  table[7] = 0.0;                          // 01 XX
}
void UpdateLookupTable(const double *af, double *lookup_table) {
  std::array<std::array<double, 8>, 5> tables;
  for (size_t i = 0; i < 4; ++i) {
    auto *lkt = lookup_table + i * 32768;
    for (size_t j = 0; j < 5; j++) {
      CalcSnpCorr(af[5 * i + j], tables[j]);
    }
    for (size_t j = 0; j < 32768; ++j) {
      lkt[j] = tables[0][j & 7u] + tables[1][j >> 3 & 7u] +
               tables[2][j >> 6 & 7u] + tables[3][j >> 9 & 7u] +
               tables[4][j >> 12 & 7u];
    }
  }
}
void UpdateLookupTable(const double *af, size_t num_snps_left,
                       double *lookup_table) {
  std::array<std::array<double, 8>, 5> tables;
  size_t shorts_left = num_snps_left / 5;
  size_t snps_remain = num_snps_left % 5;
  for (size_t i = 0; i < shorts_left; ++i) {
    auto *lkt = lookup_table + i * 32768;
    for (size_t j = 0; j < 5; j++) {
      CalcSnpCorr(af[5 * i + j], tables[j]);
    }
    for (size_t j = 0; j < 32768; ++j) {
      lkt[j] = tables[0][j & 7u] + tables[1][j >> 3 & 7u] +
               tables[2][j >> 6 & 7u] + tables[3][j >> 9 & 7u] +
               tables[4][j >> 12 & 7u];
    }
  }
  if (snps_remain != 0u) {
    auto *lkt = lookup_table + shorts_left * 32768;
    for (size_t j = 0; j < snps_remain; j++) {
      CalcSnpCorr(af[5 * shorts_left + j], tables[j]);
    }
    for (size_t j = snps_remain; j < 5; ++j) {
      std::fill(tables[j].begin(), tables[j].end(), 0.0);
    }
    for (size_t j = 0; j < 32768; ++j) {
      lkt[j] = tables[0][j & 7u] + tables[1][j >> 3 & 7u] +
               tables[2][j >> 6 & 7u] + tables[3][j >> 9 & 7u] +
               tables[4][j >> 12 & 7u];
    }
  }
}
void Mask(const uint64_t *geno64, size_t num_samples, uint64_t *mask64) {
  for (size_t i = 0; i < num_samples; ++i) {
    uint64_t tmp_geno = ~(geno64[i] ^ kMask3);
    tmp_geno = (tmp_geno & (tmp_geno >> 1)) & kMask3;
    mask64[i] = tmp_geno * 7ull;
  }
}
void UpdateMatrix(const uint64_t *geno64, const uint64_t *mask64,
                  const double *lookup_table, size_t num_samples,
                  double *matrix) {
  for (size_t i = 0; i < num_samples; i++) {
    if (mask64[i] == 0u) {
      for (size_t j = i; j < num_samples; ++j) {
        uint64_t g = (geno64[i] + geno64[j]) | mask64[j];
        matrix[i * num_samples + j] +=
            lookup_table[g & 0x7fffu] +
            lookup_table[32768 + (g >> 16 & 0x7fffu)] +
            lookup_table[65536 + (g >> 32 & 0x7fffu)] +
            lookup_table[98304 + (g >> 48 & 0x7fffu)];
      }
    } else {
      for (size_t j = i; j < num_samples; ++j) {
        uint64_t g = (geno64[i] + geno64[j]) | (mask64[i] | mask64[j]);
        matrix[i * num_samples + j] +=
            lookup_table[g & 0x7fffu] +
            lookup_table[32768 + (g >> 16 & 0x7fffu)] +
            lookup_table[65536 + (g >> 32 & 0x7fffu)] +
            lookup_table[98304 + (g >> 48 & 0x7fffu)];
      }
    }
  }
}
}  // namespace grm

namespace snplib {
void CalcIBSMatrixThread(uint8_t *geno, size_t num_samples, size_t num_snps,
                         uint64_t *matrix) {
  SNP snp(geno, num_samples);
  auto num_blocks = num_snps / 1024u;
  auto num_snps_left = num_snps % 1024u;
  auto *geno64 = new uint64_t[32 * num_samples];
  auto *mask64 = new uint64_t[32 * num_samples];
  std::fill(matrix, matrix + num_samples * num_samples, 0ull);
  for (size_t i = 0; i < num_blocks; ++i) {
    for (size_t j = 0; j < 32; ++j) {
      snp.TransposeUGeno(32, j, geno64);
      snp += 32;
    }
    ibs::Mask(geno64, num_samples, mask64);
    ibs::UpdateMatrix(geno64, mask64, num_samples, matrix);
  }
  if (num_snps_left > 0) {
    num_blocks = num_snps_left / 32;
    std::fill(geno64, geno64 + 32 * num_samples, kMask1);
    for (size_t j = 0; j < num_blocks; ++j) {
      snp.TransposeUGeno(32, j, geno64);
      snp += 32;
    }
    num_snps_left %= 32;
    if (num_snps_left > 0u) {
      snp.TransposeUGeno(num_snps_left, num_blocks, geno64);
    }
    ibs::Mask(geno64, num_samples, mask64);
    ibs::UpdateMatrix(geno64, mask64, num_samples, matrix);
  }
  delete[] geno64;
  delete[] mask64;
}
void CalcIBSConnectThread(uint8_t *src_geno, size_t num_src_samples,
                          uint8_t *dest_geno, size_t num_dest_samples,
                          size_t num_snps, uint64_t *connection) {
  SNP src_snp(src_geno, num_src_samples);
  SNP dest_snp(dest_geno, num_dest_samples);
  auto num_blocks = num_snps / 1024u;
  auto num_snps_left = num_snps % 1024u;
  auto *src_geno64 = new uint64_t[32 * num_src_samples];
  auto *src_mask64 = new uint64_t[32 * num_src_samples];
  auto *dest_geno64 = new uint64_t[32 * num_dest_samples];
  auto *dest_mask64 = new uint64_t[32 * num_dest_samples];
  std::fill(connection, connection + num_dest_samples, 0ull);
  for (size_t i = 0; i < num_blocks; ++i) {
    for (size_t j = 0; j < 32; ++j) {
      src_snp.TransposeUGeno(32, j, src_geno64);
      dest_snp.TransposeUGeno(32, j, dest_geno64);
      src_snp += 32;
      dest_snp += 32;
    }
    ibs::Mask(src_geno64, num_src_samples, src_mask64);
    ibs::Mask(dest_geno64, num_dest_samples, dest_mask64);
    ibs::UpdateConnect(src_geno64, src_mask64, dest_geno64, dest_mask64,
                       num_src_samples, num_dest_samples, connection);
  }
  if (num_snps_left > 0u) {
    std::fill(src_geno64, src_geno64 + 32 * num_src_samples, kMask1);
    std::fill(dest_geno64, dest_geno64 + 32 * num_dest_samples, kMask1);
    num_blocks = num_snps_left / 32;
    for (size_t j = 0; j < num_blocks; ++j) {
      src_snp.TransposeUGeno(32, j, src_geno64);
      dest_snp.TransposeUGeno(32, j, dest_geno64);
      src_snp += 32;
      dest_snp += 32;
    }
    num_snps_left %= 32;
    if (num_snps_left > 0u) {
      src_snp.TransposeUGeno(num_snps_left, num_blocks, src_geno64);
      dest_snp.TransposeUGeno(num_snps_left, num_blocks, dest_geno64);
    }
    ibs::Mask(src_geno64, num_src_samples, src_mask64);
    ibs::Mask(dest_geno64, num_dest_samples, dest_mask64);
    ibs::UpdateConnect(src_geno64, src_mask64, dest_geno64, dest_mask64,
                       num_src_samples, num_dest_samples, connection);
  }
  delete[] src_geno64;
  delete[] src_mask64;
  delete[] dest_geno64;
  delete[] dest_mask64;
}
void CalcGRMMatrixThread(uint8_t *geno, const double *af, size_t num_samples,
                         size_t num_snps, double *matrix) {
  SNP snp(geno, num_samples);
  auto num_blocks = num_snps / 20u;
  auto num_snps_left = num_snps % 20u;
  auto *geno64 = new uint64_t[num_samples];
  auto *mask64 = new uint64_t[num_samples];
  auto *lookup_table = new double[4 * 32768];
  for (size_t i = 0; i < num_blocks; ++i) {
    snp.TransposeSGeno(20, geno64);
    grm::UpdateLookupTable(af, lookup_table);
    grm::Mask(geno64, num_samples, mask64);
    grm::UpdateMatrix(geno64, mask64, lookup_table, num_samples, matrix);
    snp += 20;
    af += 20;
  }
  if (num_snps_left > 0) {
    snp.TransposeSGeno(num_snps_left, geno64);
    grm::UpdateLookupTable(af, num_snps_left, lookup_table);
    grm::Mask(geno64, num_samples, mask64);
    grm::UpdateMatrix(geno64, mask64, lookup_table, num_samples, matrix);
  }
  delete[] geno64;
  delete[] mask64;
  delete[] lookup_table;
}
void CalcGCTAThread(uint8_t *geno, const double *af, size_t num_samples,
                    size_t num_snps, double *diagonal) {
  SNP snp(geno, num_samples);
  auto *geno_d = new double[num_samples];
  std::fill(geno_d, geno_d + num_samples, 0.0);
  std::fill(diagonal, diagonal + num_samples, 0.0);
  std::array<double, 4> geno_table;
  for (size_t i = 0; i < num_snps; ++i) {
    auto var = 2.0 * af[i] * (1.0 - af[i]);
    geno_table[0] = 2.0 * af[i] / var;
    geno_table[1] = 0.0;
    geno_table[2] = 0.0;
    geno_table[3] = (2.0 - 2.0 * af[i]) / var;
    snp.UnpackGeno(geno_d, geno_table);
    for (size_t j = 0; j < num_samples; ++j) {
      diagonal[j] += geno_d[j];
    }
    ++snp;
  }
  delete[] geno_d;
}
void CalcUGRMMatrixThread(uint8_t *geno, size_t num_samples, size_t num_snps,
                          int64_t *matrix) {
  SNP snp(geno, num_samples);
  auto num_blocks = num_snps / 1024u;
  auto num_snps_left = num_snps % 1024u;
  auto *geno64 = new uint64_t[32 * num_samples];
  auto *mask64 = new uint64_t[32 * num_samples];
  std::fill(matrix, matrix + num_samples * num_samples, 0ull);
  for (size_t i = 0; i < num_blocks; ++i) {
    for (size_t j = 0; j < 32; ++j) {
      snp.TransposeUGeno(32, j, geno64);
      snp += 32;
    }
    ugrm::Mask(geno64, num_samples, mask64);
    ugrm::UpdateMatrix(geno64, mask64, num_samples, matrix);
  }
  if (num_snps_left > 0) {
    num_blocks = num_snps_left / 32;
    std::fill(geno64, geno64 + 32 * num_samples, kMask1);
    for (size_t j = 0; j < num_blocks; ++j) {
      snp.TransposeUGeno(32, j, geno64);
      snp += 32;
    }
    num_snps_left %= 32;
    if (num_snps_left > 0u) {
      snp.TransposeUGeno(num_snps_left, num_blocks, geno64);
    }
    ugrm::Mask(geno64, num_samples, mask64);
    ugrm::UpdateMatrix(geno64, mask64, num_samples, matrix);
  }
  delete[] geno64;
  delete[] mask64;
}
void CalcKINGThread(uint8_t *geno, size_t num_samples, size_t num_snps,
                    int64_t *matrix, uint64_t *vector) {
  SNP snp(geno, num_samples);
  auto num_blocks = num_snps / 1024u;
  auto num_snps_left = num_snps % 1024u;
  auto *geno64 = new uint64_t[32 * num_samples];
  std::fill(matrix, matrix + num_samples * num_samples, 0ull);
  std::fill(vector, vector + num_samples, 0ull);
  for (size_t i = 0; i < num_blocks; ++i) {
    for (size_t j = 0; j < 32; ++j) {
      snp.TransposeUGeno(32, j, geno64);
      snp += 32;
    }
    king::UpdateHetero(geno64, num_samples, vector);
    king::UpdateMatrix(geno64, num_samples, matrix);
  }
  if (num_snps_left > 0) {
    num_blocks = num_snps_left / 32;
    std::fill(geno64, geno64 + 32 * num_samples, kMask1);
    for (size_t j = 0; j < num_blocks; ++j) {
      snp.TransposeUGeno(32, j, geno64);
      snp += 32;
    }
    num_snps_left %= 32;
    if (num_snps_left > 0u) {
      snp.TransposeUGeno(num_snps_left, num_blocks, geno64);
    }
    king::UpdateHetero(geno64, num_samples, vector);
    king::UpdateMatrix(geno64, num_samples, matrix);
  }
  delete[] geno64;
}
void CalcAdjustedGRMThread(uint8_t *geno, size_t num_samples, size_t num_snps,
                           double *covariates, size_t num_covariates,
                           double *matrix, double *gcta_diag) {
  SNP snp(geno, num_samples);
  auto *geno_d = new double[num_samples];
  auto *mask_d = new double[num_samples];
  snplib::LogisticRegress<2> worker(num_samples, num_covariates);
  auto m = static_cast<int32_t>(num_samples);
  std::fill(matrix, matrix + num_samples * num_samples, 0.0);
  std::fill(gcta_diag, gcta_diag + num_samples, 0.0);
  for (size_t l = 0; l < num_snps; ++l) {
    snp.UnpackGeno(geno_d, mask_d);
    worker.Estimate(covariates, geno_d, mask_d);
    auto *u = worker.GetU();
    auto *w = worker.GetW();
    for (size_t i = 0; i < num_samples; ++i) {
      geno_d[i] -= u[i];
      geno_d[i] /= w[i];
      geno_d[i] *= mask_d[i];
    }
    cblas_dsyr(CblasColMajor, CblasLower, m, 1.0, geno_d, 1, matrix, m);
    for (size_t i = 0; i < num_samples; ++i) {
      gcta_diag[i] += (u[i] - 1.0) * geno_d[i] / w[i];
    }
    ++snp;
  }
  delete[] geno_d;
  delete[] mask_d;
}
void CalcAdmixedGRMthread(uint8_t *geno, size_t num_samples, size_t num_snps,
                          double *pop_af, double *pop, size_t num_pops,
                          double *matrix, double *gcta_diag) {
  SNP snp(geno, num_samples);
  auto *geno_d = new double[num_samples];
  auto *mask_d = new double[num_samples];
  auto *u = new double[num_samples];
  auto *w = new double[num_samples];
  auto d = static_cast<int32_t>(num_pops);
  auto m = static_cast<int32_t>(num_samples);
  std::fill(matrix, matrix + num_samples * num_samples, 0.0);
  std::fill(gcta_diag, gcta_diag + num_samples, 0.0);
  for (size_t l = 0; l < num_snps; ++l) {
    snp.UnpackGeno(geno_d, mask_d);
    auto *af = pop_af + l * num_pops;
    cblas_dgemv(CblasColMajor, CblasNoTrans, m, d, 2.0, pop, m, af, 1, 0.0, u,
                1);
    for (size_t i = 0; i < num_samples; ++i) {
      w[i] = std::sqrt(u[i] * (1.0 - u[i] / 2.0));
      geno_d[i] -= u[i];
      geno_d[i] /= w[i];
      geno_d[i] *= mask_d[i];
    }
    cblas_dsyr(CblasColMajor, CblasLower, m, 1.0, geno_d, 1, matrix, m);
    for (size_t i = 0; i < num_samples; ++i) {
      gcta_diag[i] += (u[i] - 1.0) * geno_d[i] / w[i];
    }
    ++snp;
  }
  delete[] geno_d;
  delete[] mask_d;
  delete[] u;
  delete[] w;
}
}  // namespace snplib