#include "grm.h"

namespace {
const uint64_t kMask = 0x1249124912491249ull;  // 0x1249=b0001001001001001
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
    uint64_t tmp_geno = ~(geno64[i] ^ kMask);
    tmp_geno = (tmp_geno & (tmp_geno >> 1)) & kMask;
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
void CalcGRMMatrixThread(uint8_t *geno, const double *af, size_t num_samples,
                         size_t num_snps, double *matrix) {
  snplib::SNP snp(geno, num_samples);
  auto num_blocks = num_snps / 20u;
  auto num_snps_left = num_snps % 20u;
  auto *geno64 = new uint64_t[num_samples];
  auto *mask64 = new uint64_t[num_samples];
  auto *lookup_table = new double[4 * 32768];
  for (size_t i = 0; i < num_blocks; ++i) {
    snplib::TransposeSGeno(snp, 20, geno64);
    UpdateLookupTable(af, lookup_table);
    Mask(geno64, num_samples, mask64);
    UpdateMatrix(geno64, mask64, lookup_table, num_samples, matrix);
    snp += 20;
    af += 20;
  }
  if (num_snps_left > 0) {
    snplib::TransposeSGeno(snp, num_snps_left, geno64);
    UpdateLookupTable(af, num_snps_left, lookup_table);
    Mask(geno64, num_samples, mask64);
    UpdateMatrix(geno64, mask64, lookup_table, num_samples, matrix);
  }
  delete[] geno64;
  delete[] mask64;
  delete[] lookup_table;
}
void CalcGCTAThread(uint8_t *geno, const double *af, size_t num_samples,
                    size_t num_snps, double *diagonal) {
  snplib::SNP snp(geno, num_samples);
  auto *geno_d = new double[num_samples];
  std::fill(geno_d, geno_d + num_samples, 0.0);
  std::fill(diagonal, diagonal + num_samples, 0.0);
  for (size_t i = 0; i < num_snps; ++i) {
    snp.UnpackGCTA(geno_d, af[i]);
    for (size_t j = 0; j < num_samples; ++j) {
      diagonal[j] += geno_d[j];
    }
    ++snp;
  }
  delete[] geno_d;
}
}  // namespace

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
void CalcGCTADiagonal(uint8_t *geno, const double *af, size_t num_samples,
                      size_t num_snps, double *diagonal, size_t num_threads) {
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