#include "statistics.h"

namespace {
std::atomic_size_t ind;
const uint64_t kMask1 = 0x5555555555555555ull;  // 0x5555=b0101010101010101
const uint64_t kMask2 = 0xAAAAAAAAAAAAAAAAull;  // 0xAAAA=b1010101010101010
const std::array<double, 4> geno_table{0.0, 0.0, 1.0, 2.0};
const std::array<double, 4> mask_table{1.0, 0.0, 1.0, 1.0};

double CalcLDR2(const uint64_t *geno_1, const uint64_t *geno_2, double af1,
                double af2, size_t num_words) {
  uint64_t num_samples = 0;
  uint64_t n11 = 0;
  uint64_t n12 = 0;
  uint64_t n21 = 0;
  uint64_t n22 = 0;
  uint64_t n_homo_1 = 0;
  uint64_t n_homo_2 = 0;
  uint64_t n_hetero_1 = 0;
  uint64_t n_hetero_2 = 0;
  for (size_t i = 0; i < num_words; ++i) {
    uint64_t mask_1 = geno_1[i] ^ kMask1;
    mask_1 = (mask_1 | (mask_1 >> 1)) & kMask1;
    uint64_t mask_2 = geno_2[i] ^ kMask1;
    mask_2 = (mask_2 | (mask_2 >> 1)) & kMask1;
    auto mask = mask_1 & mask_2;
    mask *= 3ull;
    uint64_t homo_1 = ((geno_1[i] ^ kMask2) & kMask1) & mask;
    uint64_t homo_2 = ((geno_2[i] ^ kMask2) & kMask1) & mask;
    uint64_t hetero_1 = geno_1[i] ^ kMask1;
    hetero_1 = ((hetero_1 & (hetero_1 >> 1)) & kMask1) & mask;
    uint64_t hetero_2 = geno_2[i] ^ kMask1;
    hetero_2 = ((hetero_2 & (hetero_2 >> 1)) & kMask1) & mask;
    num_samples += _mm_popcnt_u64(mask);
    mask = homo_1 & homo_2;
    n11 += _mm_popcnt_u64(mask);
    mask = homo_1 & hetero_2;
    n12 += _mm_popcnt_u64(mask);
    mask = hetero_1 & homo_2;
    n21 += _mm_popcnt_u64(mask);
    mask = hetero_1 & hetero_2;
    n22 += _mm_popcnt_u64(mask);
    n_homo_1 += _mm_popcnt_u64(homo_1);
    n_homo_2 += _mm_popcnt_u64(homo_2);
    n_hetero_1 += _mm_popcnt_u64(hetero_1);
    n_hetero_2 += _mm_popcnt_u64(hetero_2);
  }
  num_samples /= 2;
  auto n13 = n_homo_1 - n11 - n12;
  auto n23 = n_hetero_1 - n21 - n22;
  auto n31 = n_homo_2 - n11 - n21;
  auto n32 = n_hetero_2 - n12 - n22;
  auto n33 = num_samples - n11 - n12 - n13 - n21 - n22 - n23 - n31 - n32;
  auto r2 = static_cast<double>(n11) * (2.0 - 2.0 * af1) * (2.0 - 2.0 * af2);
  r2 += static_cast<double>(n21) * (1.0 - 2.0 * af1) * (2.0 - 2.0 * af2);
  r2 -= static_cast<double>(n31) * 2.0 * af1 * (2.0 - 2.0 * af2);
  r2 += static_cast<double>(n12) * (2.0 - 2.0 * af1) * (1.0 - 2.0 * af2);
  r2 += static_cast<double>(n22) * (1.0 - 2.0 * af1) * (1.0 - 2.0 * af2);
  r2 -= static_cast<double>(n32) * 2.0 * af1 * (1.0 - 2.0 * af2);
  r2 -= static_cast<double>(n13) * (2.0 - 2.0 * af1) * 2.0 * af2;
  r2 -= static_cast<double>(n23) * (1.0 - 2.0 * af1) * 2.0 * af2;
  r2 += static_cast<double>(n32) * 2.0 * af1 * 2.0 * af2;
  r2 /= 2.0 * std::sqrt(af1 * af2 * (1.0 - af1) * (1.0 - af2));
  r2 /= (num_samples - 1);
  return r2 * r2;
}
}  // namespace

namespace snplib {
void CalcAlleleFrequencies(const uint8_t *geno, size_t num_samples,
                           size_t num_snps, double *af) {
  auto num_words_snp = num_samples / 32 + ((num_samples % 32) != 0u ? 1 : 0);
  auto num_bytes = num_samples / 4 + ((num_samples % 4) != 0 ? 1 : 0);
  auto correct = 4 * num_bytes - num_samples;
  correct *= 2;
  auto *geno64 = new uint64_t[num_words_snp];
  std::fill(geno64, geno64 + num_words_snp, kMask1);
  for (size_t i = 0; i < num_snps; ++i) {
    memcpy((void *)geno64, (void *)(geno + i * num_bytes),
           sizeof(uint8_t) * num_bytes);
    uint64_t alleles = 0;
    uint64_t nonmissings = 0;
    for (size_t j = 0; j < num_words_snp; ++j) {
      uint64_t tmp_g = geno64[j] ^ kMask1;
      uint64_t tmp_m = (tmp_g | (tmp_g >> 1)) & kMask1;
      tmp_m *= 3;
      nonmissings += _mm_popcnt_u64(tmp_m);
      tmp_g = tmp_m & geno64[j];
      alleles += _mm_popcnt_u64(tmp_g);
    }
    nonmissings -= correct;
    af[i] = static_cast<double>(alleles) / nonmissings;
  }
  delete[] geno64;
}
void CalcMissing(const uint8_t *geno, size_t num_samples, size_t num_snps,
                 double *missing) {
  auto num_words_snp = num_samples / 32 + ((num_samples % 32) != 0u ? 1 : 0);
  auto num_bytes = num_samples / 4 + ((num_samples % 4) != 0 ? 1 : 0);
  auto correct = 4 * num_bytes - num_samples;
  correct *= 2;
  auto *geno64 = new uint64_t[num_words_snp];
  std::fill(geno64, geno64 + num_words_snp, kMask1);
  for (size_t i = 0; i < num_snps; ++i) {
    memcpy((void *)geno64, (void *)(geno + i * num_bytes),
           sizeof(uint8_t) * num_bytes);
    uint64_t alleles = 0;
    uint64_t nonmissings = 0;
    for (size_t j = 0; j < num_words_snp; ++j) {
      uint64_t tmp_g = geno64[j] ^ kMask1;
      uint64_t tmp_m = (tmp_g | (tmp_g >> 1)) & kMask1;
      tmp_m *= 3;
      nonmissings += _mm_popcnt_u64(tmp_m);
      tmp_g = tmp_m & geno64[j];
      alleles += _mm_popcnt_u64(tmp_g);
    }
    nonmissings -= correct;
    missing[i] = 1.0 - static_cast<double>(nonmissings) / num_samples;
  }
  delete[] geno64;
}
void CalcLDscoresThread(const uint8_t *geno, const int32_t *bp,
                        const double *af, size_t num_snps, size_t num_samples,
                        size_t window_size, double r2_threshold, double *ldcv) {
  auto local_ind = ind++;
  auto ws = static_cast<int32_t>(window_size / 2);
  auto num_words = num_samples / 32 + (num_samples % 32 > 0 ? 1 : 0);
  auto *geno_1 = new uint64_t[num_words];
  auto *geno_2 = new uint64_t[num_words];
  while (local_ind < num_snps) {
    SNP snp1(geno, num_samples);
    snp1 += local_ind;
    snp1.Copy(geno_1);
    auto lower = bp[local_ind] - ws;
    lower = lower > 0 ? lower : 0;
    auto upper = bp[local_ind] + ws;
    size_t first_ind = 0;
    while (bp[first_ind] < lower) {
      first_ind++;
    }
    size_t last_ind = first_ind;
    while ((bp[last_ind] <= upper) && (last_ind <= num_snps)) {
      last_ind++;
    }
    double cv = 0.0;
    for (size_t i = first_ind; i < last_ind; ++i) {
      if (i != local_ind) {
        SNP snp2(geno, num_samples);
        snp2 += i;
        snp2.Copy(geno_2);
        auto r2 = CalcLDR2(geno_1, geno_2, af[local_ind], af[i], num_words);
        cv += r2 > r2_threshold ? r2 : 0;
      }
    }
    ldcv[local_ind] = cv;
    local_ind = ind++;
  }
  delete[] geno_1;
  delete[] geno_2;
}
void CalcLDscores(const uint8_t *geno, const int32_t *bp, const double *af,
                  size_t num_snps, size_t num_samples, size_t window_size,
                  double r2_threshold, double *ldcv, size_t num_threads) {
  std::vector<std::thread> workers;
  ind = 0;
  set_num_threads(1);
  for (size_t i = 0; i < num_threads; ++i) {
    workers.emplace_back(CalcLDscoresThread, geno, bp, af, num_snps,
                         num_samples, window_size, r2_threshold, ldcv);
  }
  for (auto &&iter : workers) {
    iter.join();
  }
}
void CalcAdjustedAFThread(const uint8_t *geno, size_t num_samples,
                          size_t num_snps, double *covariates,
                          size_t num_covariates, double *af) {
  auto *mask_d = new double[num_samples];
  auto *geno_d = new double[num_samples];
  auto local_ind = ind++;
  LogisticRegress<2> worker(num_samples, num_covariates);
  while (local_ind < num_snps) {
    auto *tmp_af = af + local_ind * num_samples;
    SNP snp(geno, num_samples);
    snp += local_ind;
    snp.UnpackGeno(geno_table, mask_table, geno_d, mask_d);
    worker.Estimate(covariates, geno_d, mask_d);
    auto *u = worker.GetU();
    for (size_t i = 0; i < num_samples; ++i) {
      tmp_af[i] = u[i] / 2.0;
    }
    local_ind = ind++;
  }
  delete[] geno_d;
  delete[] mask_d;
}
void CalcAdjustedAF(const uint8_t *geno, size_t num_samples, size_t num_snps,
                    double *covariates, size_t num_covariates, double *af,
                    size_t num_threads) {
  std::vector<std::thread> workers;
  ind = 0;
  set_num_threads(1);
  for (size_t i = 0; i < num_threads; ++i) {
    workers.emplace_back(CalcAdjustedAFThread, geno, num_samples, num_snps,
                         covariates, num_covariates, af);
  }
  for (auto &&iter : workers) {
    iter.join();
  }
}
void CalcAdjustedMAFThread(const uint8_t *geno, size_t num_samples,
                           size_t num_snps, double *covariates,
                           size_t num_covariates, double *min_maf) {
  auto *mask_d = new double[num_samples];
  auto *geno_d = new double[num_samples];
  auto local_ind = ind++;
  LogisticRegress<2> worker(num_samples, num_covariates);
  while (local_ind < num_snps) {
    SNP snp(geno, num_samples);
    snp += local_ind;
    snp.UnpackGeno(geno_table, mask_table, geno_d, mask_d);
    worker.Estimate(covariates, geno_d, mask_d);
    auto *u = worker.GetU();
    double mm = 1.0;
    for (size_t i = 0; i < num_samples; ++i) {
      auto tmp_maf = u[i] > 1.0 ? 1.0 - u[i] / 2.0 : u[i] / 2.0;
      if (tmp_maf < mm) {
        mm = tmp_maf;
      }
    }
    min_maf[local_ind] = mm;
    local_ind = ind++;
  }
  delete[] geno_d;
  delete[] mask_d;
}
void CalcAdjustedMAF(const uint8_t *geno, size_t num_samples, size_t num_snps,
                     double *covariates, size_t num_covariates, double *min_maf,
                     size_t num_threads) {
  std::vector<std::thread> workers;
  ind = 0;
  set_num_threads(1);
  for (size_t i = 0; i < num_threads; ++i) {
    workers.emplace_back(CalcAdjustedMAFThread, geno, num_samples, num_snps,
                         covariates, num_covariates, min_maf);
  }
  for (auto &&iter : workers) {
    iter.join();
  }
}
}  // namespace snplib