#include "basic_statistics.h"
namespace {
const uint64_t kMask = 0x5555555555555555;  // 0x5555=b0101010101010101
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
}  // namespace snplib