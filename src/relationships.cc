#include "relationships.h"

namespace {
const uint64_t kMask1 = 0x1249124912491249ull;  // 0x1249=b0001001001001001
const uint64_t kMask2 = 0x5555555555555555ull;  // 0x5555=b0101010101010101
}  // namespace

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
    uint64_t tmp_geno = ~(geno64[i] ^ kMask1);
    tmp_geno = (tmp_geno & (tmp_geno >> 1)) & kMask1;
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

namespace ibs {
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
    uint64_t tmp_geno = geno64[i] ^ kMask2;
    tmp_geno = (tmp_geno | (tmp_geno >> 1)) & kMask2;
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
}  // namespace ibs

namespace snplib {
void CalcAdjustedGRMThread(uint8_t *geno, size_t num_samples, size_t num_snps,
                           double *covariates, size_t num_covariates,
                           double *matrix, double *gcta_diag) {
  SNP snp(geno, num_samples);
  auto *geno_d = new double[num_samples];
  auto *mask_d = new double[num_samples];
  LogisticRegress<2> worker(num_samples, num_covariates);
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
void CalcAdjustedGRM(uint8_t *geno, size_t num_samples, size_t num_snps,
                     double *covariates, size_t num_covariates, double *matrix,
                     double *gcta_diag, size_t num_threads) {
  std::vector<std::thread> workers(num_threads);
  auto *matrices = new double[num_samples * num_samples * num_threads];
  auto *diags = new double[num_samples * num_threads];
  auto num_snps_job = num_snps / num_threads + 1;
  auto num_snps_left = num_snps % num_threads;
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + (num_samples_left > 0 ? 1 : 0);
  for (size_t i = 0; i < num_snps_left; ++i) {
    workers[i] = std::thread(CalcAdjustedGRMThread, geno, num_samples,
                             num_snps_job, covariates, num_covariates,
                             matrices + i * num_samples * num_samples,
                             diags + i * num_samples);
    geno += num_snps_job * num_bytes;
  }
  --num_snps_job;
  for (size_t i = num_snps_left; i < num_threads; ++i) {
    workers[i] = std::thread(CalcAdjustedGRMThread, geno, num_samples,
                             num_snps_job, covariates, num_covariates,
                             matrices + i * num_samples * num_samples,
                             diags + i * num_samples);
    geno += num_snps_job * num_bytes;
  }
  for (auto &&iter : workers) {
    iter.join();
  }
  std::fill(matrix, matrix + num_samples * num_samples, 0.0);
  std::fill(gcta_diag, gcta_diag + num_samples, 0.0);
  for (size_t k = 0; k < num_threads; ++k) {
    auto *tmp_m = matrices + k * num_samples * num_samples;
    auto *tmp_d = diags + k * num_samples;
    for (size_t i = 0; i < num_samples; ++i) {
      for (size_t j = i; j < num_samples; ++j) {
        matrix[i * num_samples + j] += tmp_m[i * num_samples + j];
      }
      gcta_diag[i] += tmp_d[i];
    }
  }
  delete[] matrices;
  for (size_t i = 0; i < num_samples; ++i) {
    for (size_t j = i; j < num_samples; ++j) {
      matrix[i * num_samples + j] /= num_snps;
      matrix[j * num_samples + i] = matrix[i * num_samples + j];
    }
    gcta_diag[i] /= num_snps;
  }
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
void CalcAdmixedGRM(uint8_t *geno, size_t num_samples, size_t num_snps,
                    double *pop_af, double *pop, size_t num_pops,
                    double *matrix, double *gcta_diag, size_t num_threads) {
  std::vector<std::thread> workers(num_threads);
  auto *matrices = new double[num_samples * num_samples * num_threads];
  auto *diags = new double[num_samples * num_threads];
  auto num_snps_job = num_snps / num_threads + 1;
  auto num_snps_left = num_snps % num_threads;
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + (num_samples_left > 0 ? 1 : 0);
  for (size_t i = 0; i < num_snps_left; ++i) {
    workers[i] = std::thread(CalcAdmixedGRMthread, geno, num_samples,
                             num_snps_job, pop_af, pop, num_pops,
                             matrices + i * num_samples * num_samples,
                             diags + i * num_samples);
    pop_af += num_snps_job * num_pops;
    geno += num_snps_job * num_bytes;
  }
  --num_snps_job;
  for (size_t i = num_snps_left; i < num_threads; ++i) {
    workers[i] = std::thread(CalcAdmixedGRMthread, geno, num_samples,
                             num_snps_job, pop_af, pop, num_pops,
                             matrices + i * num_samples * num_samples,
                             diags + i * num_samples);
    pop_af += num_snps_job * num_pops;
    geno += num_snps_job * num_bytes;
  }
  for (auto &&iter : workers) {
    iter.join();
  }
  std::fill(matrix, matrix + num_samples * num_samples, 0.0);
  std::fill(gcta_diag, gcta_diag + num_samples, 0.0);
  for (size_t k = 0; k < num_threads; ++k) {
    auto *tmp_m = matrices + k * num_samples * num_samples;
    auto *tmp_d = diags + k * num_samples;
    for (size_t i = 0; i < num_samples; ++i) {
      for (size_t j = i; j < num_samples; ++j) {
        matrix[i * num_samples + j] += tmp_m[i * num_samples + j];
      }
      gcta_diag[i] += tmp_d[i];
    }
  }
  delete[] matrices;
  for (size_t i = 0; i < num_samples; ++i) {
    for (size_t j = i; j < num_samples; ++j) {
      matrix[i * num_samples + j] /= num_snps;
      matrix[j * num_samples + i] = matrix[i * num_samples + j];
    }
    gcta_diag[i] /= num_snps;
  }
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
void CalcGRMMatrix(uint8_t *geno, const double *af, size_t num_samples,
                   size_t num_snps, double *grm, size_t num_threads) {
  std::vector<std::thread> workers(num_threads);
  auto *matrices = new double[num_samples * num_samples * num_threads];
  auto num_snps_job = num_snps / num_threads + 1;
  auto num_snps_left = num_snps % num_threads;
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + (num_samples_left > 0 ? 1 : 0);
  for (size_t i = 0; i < num_snps_left; ++i) {
    workers[i] =
        std::thread(CalcGRMMatrixThread, geno, af, num_samples, num_snps_job,
                    matrices + i * num_samples * num_samples);
    geno += num_snps_job * num_bytes;
    af += num_snps_job;
  }
  --num_snps_job;
  for (size_t i = num_snps_left; i < num_threads; ++i) {
    workers[i] =
        std::thread(CalcGRMMatrixThread, geno, af, num_samples, num_snps_job,
                    matrices + i * num_samples * num_samples);
    geno += num_snps_job * num_bytes;
    af += num_snps_job;
  }
  for (auto &&iter : workers) {
    iter.join();
  }
  std::fill(grm, grm + num_samples * num_samples, 0.0);
  for (size_t k = 0; k < num_threads; ++k) {
    auto *tmp_m = matrices + k * num_samples * num_samples;
    for (size_t i = 0; i < num_samples; ++i) {
      for (size_t j = i; j < num_samples; ++j) {
        grm[i * num_samples + j] += tmp_m[i * num_samples + j];
      }
    }
  }
  delete[] matrices;
  for (size_t i = 0; i < num_samples; ++i) {
    for (size_t j = i; j < num_samples; ++j) {
      grm[i * num_samples + j] /= num_snps;
      grm[j * num_samples + i] = grm[i * num_samples + j];
    }
  }
}
void CalcGCTAThread(uint8_t *geno, const double *af, size_t num_samples,
                    size_t num_snps, double *diagonal) {
  SNP snp(geno, num_samples);
  auto *geno_d = new double[num_samples];
  std::fill(geno_d, geno_d + num_samples, 0.0);
  std::fill(diagonal, diagonal + num_samples, 0.0);
  std::array<double, 4> geno_table{0.0, 0.0, 0.0, 0.0};
  for (size_t i = 0; i < num_snps; ++i) {
    auto var = 2.0 * af[i] * (1.0 - af[i]);
    geno_table[0] = 2.0 * af[i] / var;
    geno_table[3] = (2.0 - 2.0 * af[i]) / var;
    snp.UnpackGeno(geno_d, geno_table);
    for (size_t j = 0; j < num_samples; ++j) {
      diagonal[j] += geno_d[j];
    }
    ++snp;
  }
  delete[] geno_d;
}
void CalcGCTADiagonal(uint8_t *geno, const double *af, size_t num_samples,
                      size_t num_snps, double *diagnol, size_t num_threads) {
  std::vector<std::thread> workers(num_threads);
  auto *diagonals = new double[num_samples * num_threads];
  auto num_snps_job = num_snps / num_threads + 1;
  auto num_snps_left = num_snps % num_threads;
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + (num_samples_left > 0 ? 1 : 0);
  for (size_t i = 0; i < num_snps_left; ++i) {
    workers[i] = std::thread(CalcGCTAThread, geno, af, num_samples,
                             num_snps_job, diagonals + i * num_samples);
    geno += num_snps_job * num_bytes;
    af += num_snps_job;
  }
  --num_snps_job;
  for (size_t i = num_snps_left; i < num_threads; ++i) {
    workers[i] = std::thread(CalcGCTAThread, geno, af, num_samples,
                             num_snps_job, diagonals + i * num_samples);
    geno += num_snps_job * num_bytes;
    af += num_snps_job;
  }
  for (auto &&iter : workers) {
    iter.join();
  }
  std::fill(diagonal, diagonal + num_samples, 0.0);
  for (size_t k = 0; k < num_threads; ++k) {
    auto *tmp_d = diagonals + k * num_samples;
    for (size_t i = 0; i < num_samples; ++i) {
      diagonal[i] += tmp_d[i];
    }
  }
  delete[] diagonals;
  for (size_t i = 0; i < num_samples; ++i) {
    diagonal[i] /= num_snps;
  }
}
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
    std::fill(geno64, geno64 + 32 * num_samples, kMask);
    for (size_t j = 0; j < num_blocks; ++j) {
      snp.TransposeUGeno(32, j, geno64);
      snp += 32;
    }
    num_snps_left %= 32;
    if (num_snps_left > 0u) {
      snp.TransposeUGeno(num_snps_left, num_blocks, geno64);
    }
    Mask(geno64, num_samples, mask64);
    UpdateMatrix(geno64, mask64, num_samples, matrix);
  }
  delete[] geno64;
  delete[] mask64;
}
void CalcIBSMatrix(uint8_t *geno, size_t num_samples, size_t num_snps,
                   double *ibs, size_t num_threads) {}
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
void CalcIBSConnection(uint8_t *src_geno, size_t num_src_samples,
                       uint8_t *dest_geno, size_t num_dest_samples,
                       size_t num_snps, double *connection,
                       size_t num_threads) {}
void CalcKINGMatrix(uint8_t *geno, size_t num_samples, size_t num_snps,
                    double *matrix, size_t num_threads) {}
void CalcUGRMMatrix(uint8_t *geno, size_t num_samples, size_t num_snps,
                    double *matrix, size_t num_threads) {}
}  // namespace snplib