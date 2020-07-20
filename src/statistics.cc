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
  for (size_t i = 0; i < num_words; ++i) {
    uint64_t homo_1 = (geno_1[i] ^ kMask2) & kMask1;
    uint64_t homo_2 = (geno_2[i] ^ kMask2) & kMask1;
    uint64_t hetero_1 = geno_1[i] ^ kMask1;
    hetero_1 = (hetero_1 ^ kMask2) & kMask1;
    uint64_t hetero_2 = geno_2[i] ^ kMask1;
    hetero_2 = (hetero_2 ^ kMask2) & kMask1;
    uint64_t mask_1 = geno_1[i] ^ kMask1;
    mask_1 = (mask_1 | (mask_1 >> 1)) & kMask1;
    uint64_t mask_2 = geno_2[i] ^ kMask1;
    mask_2 = (mask_2 | (mask_2 >> 1)) & kMask1;
    auto mask = mask_1 & mask_2;
    num_samples += _mm_popcnt_u64(mask);
    auto tmp_m = homo_1 & homo_2 & mask;
    n11 += _mm_popcnt_u64(tmp_m);
    tmp_m = homo_1 & hetero_2 & mask;
    n12 += _mm_popcnt_u64(tmp_m);
    tmp_m = hetero_1 & homo_2 & mask;
    n21 += _mm_popcnt_u64(tmp_m);
    tmp_m = hetero_1 & hetero_2 & mask;
    n22 += _mm_popcnt_u64(tmp_m);
  }
  auto n1 = static_cast<double>(num_samples);
  auto n2 = static_cast<double>(2 * n11 + n12 + n21);
  auto n3 = static_cast<double>(n22);
  auto pq = 1.0 - 2.0 * af1 - 2.0 * af2;
  auto a = 4.0 * n1;
  auto b = 2.0 * n1 * pq - 2.0 * n2 - n3;
  auto c = 2.0 * n1 * af1 * af2 - n2 * pq - n3 * (1.0 - af1 - af2);
  auto d = -n2 * af1 * af2;
  auto x_n = -b / a / 3.0;
  auto delta = (b * b - 3.0 * a * c) / a / a / 9.0;
  auto h = 4.0 * a * a * std::pow(delta, 6);
  auto y_n = a * std::pow(x_n, 3) + b * x_n * x_n + c * x_n + d;
  auto dd = y_n * y_n - h;
  if (dd > 0.0) {
    auto result = x_n + std::cbrt((std::sqrt(dd) - y_n) / a / 2.0) +
                  std::cbrt((-std::sqrt(dd) - y_n) / a / 2.0);
    return result;
  }
  if (dd == 0.0) {
    auto mu = std::cbrt(y_n / a / 2.0);
  }
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