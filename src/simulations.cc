#include "simulations.h"

namespace snplib {
void UpdateAf(const double *aaf, size_t num_pops, size_t num_snps,
              size_t num_generations, size_t effective_sample_size,
              double *af) {
  std::random_device rd;
  std::mt19937_64 gen(rd());
  for (size_t i = 0; i < num_pops; ++i) {
    auto *tmp_af = af + i * num_snps;
    std::copy(aaf, aaf + num_snps, tmp_af);
    for (size_t j = 0; j < num_snps; ++j) {
      auto tmp_a = tmp_af[j];
      for (size_t k = 0; k < num_generations; ++k) {
        std::binomial_distribution<> bdist(effective_sample_size, tmp_a);
        auto N = bdist(gen);
        if (N == 0) {
          N = 1;
        }
        if (N == effective_sample_size) {
          N -= 1;
        }
        tmp_a = static_cast<double>(N) / effective_sample_size;
      }
      tmp_af[j] = tmp_a;
    }
  }
}
void GenerateIndividuals(const double *af, size_t num_samples, size_t num_snps,
                         uint8_t *geno) {
  std::random_device rd;
  std::mt19937_64 gen(rd());
  std::uniform_real_distribution<double> dis(0.0, 1.0);
  auto num_bytes = num_samples / 4 + (num_samples % 4 > 0u ? 1 : 0);
  std::fill(geno, geno + num_bytes * num_snps, 0);
  for (size_t i = 0; i < num_snps; ++i) {
    auto tmp_g = geno + i * num_bytes;
    auto tmp_af = af + i * num_samples;
    for (size_t j = 0; j < num_samples; ++j) {
      uint8_t t = (dis(gen) < tmp_af[j]) + (dis(gen) < tmp_af[j]);
      t = t > 0u ? t + 1 : 0;
      auto l = j >> 2;
      auto k = (j & 3u) << 1;
      tmp_g[l] += t << k;
    }
  }
}
void GenerateAdmixedIndividuals(const double *af, size_t num_snps,
                                size_t num_samples, uint8_t *geno) {
  std::random_device rd;
  std::mt19937_64 gen(rd());
  std::uniform_real_distribution<double> dis(0.0, 1.0);
  auto num_bytes = num_samples / 4 + (num_samples % 4 > 0u ? 1 : 0);
  std::fill(geno, geno + num_bytes * num_snps, 0);
  for (size_t i = 0; i < num_snps; ++i) {
    auto tmp_g = geno + i * num_bytes;
    for (size_t j = 0; j < num_samples; ++j) {
      uint8_t t = (dis(gen) < af[2 * i]) + (dis(gen) < af[2 * i + 1]);
      t = t > 0u ? t + 1 : 0;
      auto l = j >> 2;
      auto k = (j & 3u) << 1;
      tmp_g[l] += t << k;
    }
  }
}
void GeneratePairwiseSiblings(uint8_t *parent_geno, size_t num_families,
                              size_t num_snps, uint8_t *siblings_geno) {
  std::random_device rd;
  std::mt19937_64 gen(rd());
  std::bernoulli_distribution bdis(0.5);
  auto num_samples = num_families * 2;
  auto num_bytes = num_samples / 4 + (num_samples % 4 > 0u ? 1 : 0);
  std::fill(siblings_geno, siblings_geno + num_bytes * num_snps, 0u);
  for (size_t i = 0; i < num_snps; ++i) {
    auto *tmp_p = parent_geno + i * num_bytes;
    auto *tmp_s = siblings_geno + i * num_bytes;
    for (size_t j = 0; j < num_families; ++j) {
      auto l = j / 2;
      auto k = (j % 2) * 4;
      uint8_t p1_g = (tmp_p[l] >> k) & 3u;
      p1_g = p1_g > 0u ? p1_g - 1 : 0;
      uint8_t p2_g = (tmp_p[l] >> (k + 2)) & 3u;
      p2_g = p2_g > 0u ? p2_g - 1 : 0;
      uint8_t t = (p1_g == 1 ? bdis(gen) : p1_g / 2) +
                  (p2_g == 1 ? bdis(gen) : p2_g / 2);
      t = t > 0u ? t + 1 : 0;
      tmp_s[l] += t << k;
      t = (p1_g == 1 ? bdis(gen) : p1_g / 2) +
          (p2_g == 1 ? bdis(gen) : p2_g / 2);
      t = t > 0u ? t + 1 : 0;
      tmp_s[l] += t << (k + 2);
    }
  }
}
}  // namespace snplib