#include "statistics.h"

namespace {
const uint64_t kMask = 0x5555555555555555;  // 0x5555=b0101010101010101
void CalcAdjustedAFThread(uint8_t *geno, size_t num_samples, size_t num_snps,
                          double *covariates, size_t num_covariates,
                          double *af) {
  auto *mask_d = new double[num_samples];
  auto *geno_d = new double[num_samples];
  snplib::SNP snp(geno, num_samples);
  snplib::LogisticRegress<2> worker(num_samples, num_covariates);
  for (size_t l = 0; l < num_snps; ++l) {
    auto *tmp_af = af + l * num_samples;
    snp.UnpackGeno(geno_d, mask_d);
    worker.Estimate(covariates, geno_d, mask_d);
    auto *u = worker.GetU();
    auto *w = worker.GetW();
    double mm = 1.0;
    for (size_t i = 0; i < num_samples; ++i) {
      tmp_af[i] = u[i] / 2.0;
    }
    ++snp;
  }
  delete[] mask_d;
  delete[] geno_d;
}
void CalcAdjustedMAFThread(uint8_t *geno, size_t num_samples, size_t num_snps,
                           double *covariates, size_t num_covariates,
                           double *min_maf) {
  auto *mask_d = new double[num_samples];
  auto *geno_d = new double[num_samples];
  snplib::SNP snp(geno, num_samples);
  snplib::LogisticRegress<2> worker(num_samples, num_covariates);
  for (size_t l = 0; l < num_snps; ++l) {
    snp.UnpackGeno(geno_d, mask_d);
    worker.Estimate(covariates, geno_d, mask_d);
    auto *u = worker.GetU();
    auto *w = worker.GetW();
    double mm = 1.0;
    for (size_t i = 0; i < num_samples; ++i) {
      auto tmp_maf = u[i] > 1.0 ? 1.0 - u[i] / 2.0 : u[i] / 2.0;
      if (tmp_maf < mm) {
        mm = tmp_maf;
      }
    }
    min_maf[l] = mm;
    ++snp;
  }
  delete[] mask_d;
  delete[] geno_d;
}
}  // namespace

namespace snplib {
void CalcAlleleFrequencies(uint8_t *geno, size_t num_samples, size_t num_snps,
                           double *af) {
  SNP snp(geno, num_samples);
  auto num_words_snp = num_samples / 32 + ((num_samples % 32) != 0u ? 1 : 0);
  auto num_bytes = num_samples / 4 + ((num_samples % 4) != 0 ? 1 : 0);
  auto correct = 4 * num_bytes - num_samples;
  correct *= 2;
  auto *geno64 = new uint64_t[num_words_snp];
  std::fill(geno64, geno64 + num_words_snp, kMask);
  for (size_t i = 0; i < num_snps; ++i) {
    snp.Copy(geno64);
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
    ++snp;
  }
  delete[] geno64;
}
void CalcMissing(uint8_t *geno, size_t num_samples, size_t num_snps,
                 double *ms) {
  SNP snp(geno, num_samples);
  auto num_words_snp = num_samples / 32 + ((num_samples % 32) != 0u ? 1 : 0);
  auto num_bytes = num_samples / 4 + ((num_samples % 4) != 0 ? 1 : 0);
  auto correct = 4 * num_bytes - num_samples;
  correct *= 2;
  auto *geno64 = new uint64_t[num_words_snp];
  std::fill(geno64, geno64 + num_words_snp, kMask);
  for (size_t i = 0; i < num_snps; ++i) {
    snp.Copy(geno64);
    uint64_t nonmissings = 0;
    for (size_t j = 0; j < num_words_snp; ++j) {
      uint64_t tmp_g = geno64[j] ^ kMask;
      uint64_t tmp_m = (tmp_g | (tmp_g >> 1)) & kMask;
      tmp_m *= 3;
      nonmissings += _mm_popcnt_u64(tmp_m);
    }
    nonmissings -= correct;
    ms[i] = 1.0 - static_cast<double>(nonmissings) / num_samples / 2.0;
    ++snp;
  }
  delete[] geno64;
}
void CalcAdjustedAF(uint8_t *geno, size_t num_samples, size_t num_snps,
                    double *covariates, size_t num_covariates, double *af,
                    size_t num_threads) {
  std::vector<std::thread> workers(num_threads);
  auto num_snps_job = num_snps / num_threads + 1;
  auto num_snps_left = num_snps % num_threads;
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + (num_samples_left > 0 ? 1 : 0);
  for (size_t i = 0; i < num_snps_left; ++i) {
    workers[i] = std::thread(CalcAdjustedAFThread, geno, num_samples,
                             num_snps_job, covariates, num_covariates, af);
    geno += num_bytes * num_snps_job;
    af += num_samples * num_snps_job;
  }
  --num_snps_job;
  for (size_t i = num_snps_left; i < num_threads; ++i) {
    workers[i] = std::thread(CalcAdjustedAFThread, geno, num_samples,
                             num_snps_job, covariates, num_covariates, af);
    geno += num_bytes * num_snps_job;
    af += num_samples * num_snps_job;
  }
  for (auto &&iter : workers) {
    iter.join();
  }
}
void CalcAdjustedMAF(uint8_t *geno, size_t num_samples, size_t num_snps,
                     double *covariates, size_t num_covariates, double *min_maf,
                     size_t num_threads) {
  std::vector<std::thread> workers(num_threads);
  auto num_snps_job = num_snps / num_threads + 1;
  auto num_snps_left = num_snps % num_threads;
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + (num_samples_left > 0 ? 1 : 0);
  for (size_t i = 0; i < num_snps_left; ++i) {
    workers[i] = std::thread(CalcAdjustedMAFThread, geno, num_samples,
                             num_snps_job, covariates, num_covariates, min_maf);
    geno += num_bytes * num_snps_job;
    min_maf += num_snps_job;
  }
  --num_snps_job;
  for (size_t i = num_snps_left; i < num_threads; ++i) {
    workers[i] = std::thread(CalcAdjustedMAFThread, geno, num_samples,
                             num_snps_job, covariates, num_covariates, min_maf);
    geno += num_bytes * num_snps_job;
    min_maf += num_snps_job;
  }
  for (auto &&iter : workers) {
    iter.join();
  }
}
}  // namespace snplib