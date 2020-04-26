#include "ancestry_estimation.h"

namespace snplib {
void CalcPCALoadingsExact(uint8_t *geno, size_t num_samples, size_t num_snps,
                          size_t num_components, double *loadings,
                          size_t num_threads) {
  set_num_threads(num_threads);
  auto *af = new double[num_snps];
  auto *geno_d = new double[num_samples * num_snps];
  SNP snp(geno, num_samples);
  std::array<double, 4> geno_table{0.0, 0.0, 0.0, 0.0};
  CalcAlleleFrequencies(geno, num_samples, num_snps, af);
  for (size_t i = 0; i < num_snps; ++i) {
    auto std = std::sqrt(2.0 * af[i] * (1.0 - af[i]));
    auto mu = 2.0 * af[i];
    geno_table = {-mu / std, 0.0, (1.0 - mu) / std, (2.0 - mu) / std};
    auto *tmp_geno_d = geno_d + i * num_samples;
    snp.UnpackGeno(tmp_geno_d, geno_table);
    ++snp;
  }
  int32_t m = static_cast<int32_t>(num_samples);
  int32_t n = static_cast<int32_t>(num_snps);
  auto *u = new double[num_samples * num_samples];
  LAPACKE_dgesdd(LAPACK_COL_MAJOR, 'O', m, n, geno_d, m, af, u, m, nullptr, m);
  for (size_t i = 0; i < num_snps; ++i) {
    auto *tmp_geno_d = geno_d + i * num_samples;
    auto *tmp_loadings = loadings + i * num_components;
    for (size_t j = 0; j < num_components; ++j) {
      tmp_loadings[j] = tmp_geno_d[j] / af[j];
    }
  }
  delete[] af;
  delete[] geno_d;
  delete[] u;
}
void CalcPCALoadingsApprox(uint8_t *geno, size_t num_samples, size_t num_snps,
                           size_t num_components, double *loadings,
                           size_t num_parts, size_t num_threads) {
  std::random_device rd;
  std::mt19937_64 gen(rd());
  std::normal_distribution<double> dist;
  auto L = 2 * num_components;
  const size_t I = 10;
  auto *G = new double[num_samples * (I + 1) * L];
  for (size_t i = 0; i < num_samples * L; ++i) {
    G[i] = dist(gen);
  }
  auto m = static_cast<int32_t>(num_samples);
  auto n = static_cast<int32_t>(L);
  auto *af = new double[num_snps];
  auto *grm = new double[num_samples * num_samples];
  CalcAlleleFrequencies(geno, num_samples, num_snps, af);
  CalcGRMMatrix(geno, af, num_samples, num_snps, grm, num_threads);
  set_num_threads(num_threads);
  for (size_t i = 0; i < I; ++i) {
    auto *this_G = G + (i + 1) * num_samples * L;
    auto *prev_G = G + i * num_samples * L;
    cblas_dsymm(CblasColMajor, CblasLeft, CblasLower, m, n, 1.0, grm, m, prev_G,
                m, 0.0, this_G, m);
  }
  auto num_snps_part = num_snps / num_parts;
  auto *H = new double[num_snps * L * (I + 1)];
  auto *A = new double[num_samples * num_snps_part];
  n = static_cast<int32_t>(L * (I + 1));
  SNP snp(geno, num_samples);
  std::array<double, 4> geno_table{0.0, 0.0, 0.0, 0.0};
  for (size_t i = 0; i < num_parts; ++i) {
    auto *tmp_af = af + i * num_snps_part;
    auto *tmp_A = A + i * num_snps_part * num_samples;
    auto *tmp_H = H + i * num_snps_part;
    for (size_t j = 0; j < num_snps_part; ++j) {
      auto std = std::sqrt(2.0 * tmp_af[j] * (1.0 - tmp_af[j]));
      auto mu = 2.0 * tmp_af[j];
      geno_table = {-mu / std, 0.0, (1.0 - mu) / std, (2.0 - mu) / std};
      snp.UnpackGeno(tmp_A + j * num_samples, geno_table);
      ++snp;
    }
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, )
  }
  delete[] grm;
  delete[] af;
  delete[] H;
  delete[] A;
}
void ProjectPCA(uint8_t *dest_geno, size_t num_samples, size_t num_snps,
                double *loadings, size_t num_components, double *scores,
                size_t num_parts, size_t num_threads) {}
void CalcSUGIBSLoadingsExact(uint8_t *geno, size_t num_samples, size_t num_snps,
                             size_t num_components, double *loadings,
                             size_t num_threads) {}
void CalcSUGIBSLoadingsApprox(uint8_t *geno, size_t num_samples,
                              size_t num_snps, size_t num_components,
                              double *loadings, size_t num_parts,
                              size_t num_threads) {}
void ProjectSUGIBS(uint8_t *src_geno, size_t num_src_samples,
                   uint8_t *dest_geno, size_t num_dest_samples, size_t num_snps,
                   double *loadings, size_t num_components, double *scores,
                   size_t num_parts, size_t num_threads) {}
void CalcUPCALoadingsExact(uint8_t *geno, size_t num_samples, size_t num_snps,
                           size_t num_components, double *loadings,
                           size_t num_threads) {
  const std::array<double, 4> geno_table{-1.0, 0.0, 0.0, 1.0};
  set_num_threads(num_threads);
  auto *geno_d = new double[num_samples * num_snps];
  SNP snp(geno, num_samples);
  for (size_t i = 0; i < num_snps; ++i) {
    auto *tmp_geno_d = geno_d + i * num_samples;
    snp.UnpackGeno(tmp_geno_d, geno_table);
    ++snp;
  }
  int32_t m = static_cast<int32_t>(num_samples);
  int32_t n = static_cast<int32_t>(num_snps);
  auto *u = new double[num_samples * num_samples];
  auto *s = new double[num_samples];
  LAPACKE_dgesdd(LAPACK_COL_MAJOR, 'O', m, n, geno_d, m, s, u, m, nullptr, m);
  for (size_t i = 0; i < num_snps; ++i) {
    auto *tmp_geno_d = geno_d + i * num_samples;
    auto *tmp_loadings = loadings + i * num_components;
    for (size_t j = 0; j < num_components; ++j) {
      tmp_loadings[j] = tmp_geno_d[j] / s[j];
    }
  }
  delete[] s;
  delete[] geno_d;
  delete[] u;
}

void CalcUPCALoadingsApprox(uint8_t *geno, size_t num_samples, size_t num_snps,
                            size_t num_components, double *loadings,
                            size_t num_parts, size_t num_threads) {}
void ProjectUPCA(uint8_t *dest_geno, size_t num_samples, size_t num_snps,
                 double *loadings, size_t num_components, double *scores,
                 size_t num_parts, size_t num_threads) {}
}  // namespace snplib