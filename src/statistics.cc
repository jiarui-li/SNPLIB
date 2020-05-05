#include "statistics.h"

namespace {
std::atomic_size_t ind;
const uint64_t kMask = 0x5555555555555555;  // 0x5555=b0101010101010101
const std::array<double, 4> geno_table{0.0, 0.0, 1.0, 2.0};
const std::array<double, 4> mask_table{1.0, 0.0, 1.0, 1.0};
}  // namespace

namespace snplib {
void CalcAlleleFrequencies(const uint8_t *geno, size_t num_samples,
                           size_t num_snps, double *af) {
  auto num_words_snp = num_samples / 32 + ((num_samples % 32) != 0u ? 1 : 0);
  auto num_bytes = num_samples / 4 + ((num_samples % 4) != 0 ? 1 : 0);
  auto correct = 4 * num_bytes - num_samples;
  correct *= 2;
  auto *geno64 = new uint64_t[num_words_snp];
  std::fill(geno64, geno64 + num_words_snp, kMask);
  for (size_t i = 0; i < num_snps; ++i) {
    memcpy((void *)geno64, (void *)(geno + i * num_bytes),
           sizeof(uint8_t) * num_bytes);
    uint64_t alleles = 0;
    uint64_t nonmissings = 0;
    for (size_t j = 0; j < num_words_snp; ++j) {
      uint64_t tmp_g = geno64[j] ^ kMask;
      uint64_t tmp_m = (tmp_g | (tmp_g >> 1)) & kMask;
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
  std::fill(geno64, geno64 + num_words_snp, kMask);
  for (size_t i = 0; i < num_snps; ++i) {
    memcpy((void *)geno64, (void *)(geno + i * num_bytes),
           sizeof(uint8_t) * num_bytes);
    uint64_t alleles = 0;
    uint64_t nonmissings = 0;
    for (size_t j = 0; j < num_words_snp; ++j) {
      uint64_t tmp_g = geno64[j] ^ kMask;
      uint64_t tmp_m = (tmp_g | (tmp_g >> 1)) & kMask;
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
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + num_samples_left > 0 ? 1 : 0;
  auto *mask_d = new double[num_samples];
  auto *geno_d = new double[num_samples];
  auto local_ind = ind++;
  LogisticRegress<2> worker(num_samples, num_covariates);
  while (local_ind < num_snps) {
    auto *tmp_af = af + local_ind * num_samples;
    auto *tmp_g = geno + local_ind * num_bytes;
    UnpackGeno(tmp_g, num_samples, geno_table, mask_table, geno_d, mask_d);
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
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + num_samples_left > 0 ? 1 : 0;
  auto *mask_d = new double[num_samples];
  auto *geno_d = new double[num_samples];
  auto local_ind = ind++;
  LogisticRegress<2> worker(num_samples, num_covariates);
  while (local_ind < num_snps) {
    auto *tmp_g = geno + local_ind * num_bytes;
    UnpackGeno(tmp_g, num_samples, geno_table, mask_table, geno_d, mask_d);
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