#include "base.h"

namespace {
const uint64_t kMask = 0x5555555555555555;  // 0x5555=b0101010101010101
}  // namespace

void CalcAlleleFrequencies(uint8_t *geno, size_t num_samples, size_t num_snps,
                           double *af) {
  snplib::SNP snp(geno, num_samples);
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
void CalcCallrates(uint8_t *geno, size_t num_samples, size_t num_snps,
                   double *cr) {
  snplib::SNP snp(geno, num_samples);
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
    cr[i] = static_cast<double>(nonmissings) / num_samples / 2.0;
    ++snp;
  }
  delete[] geno64;
}
void FlipGeno(uint8_t *geno, size_t num_samples, size_t num_snps,
              const int32_t *idx) {
  snplib::SNP snp(geno, num_samples);
  for (size_t i = 0; i < num_snps; ++i) {
    auto tmp_s = snp;
    tmp_s += idx[i];
    tmp_s.Flip();
  }
}
void Keep(uint8_t *src_geno, uint8_t *dest_geno, size_t num_src_samples,
          size_t num_dest_samples, size_t num_snps, const int32_t *idx) {
  snplib::SNP src_snp(src_geno, num_src_samples);
  auto num_full_bytes = num_dest_samples / 4;
  auto num_samples_left = num_dest_samples % 4;
  auto num_bytes = num_full_bytes + (num_samples_left > 0 ? 1 : 0);
  for (size_t l = 0; l < num_snps; ++l) {
    auto *tmp_g = dest_geno + l * num_bytes;
    for (size_t i = 0; i < num_full_bytes; ++i) {
      uint8_t t = src_snp(idx[4 * i]);
      t += src_snp(idx[4 * i + 1]) << 2;
      t += src_snp(idx[4 * i + 2]) << 4;
      t += src_snp(idx[4 * i + 3]) << 6;
      tmp_g[i] = t;
    }
    if (num_samples_left > 0) {
      uint8_t t = 0;
      for (size_t i = 0; i < num_samples_left; ++i) {
        t += src_snp(idx[4 * num_full_bytes + i]) << (2 * i);
      }
      tmp_g[num_full_bytes] = t;
    }
    ++src_snp;
  }
}