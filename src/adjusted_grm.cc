#include "adjusted_grm.h"

namespace {
const std::array<double, 4> geno_table{0.0, 0.0, 1.0, 2.0};
const std::array<double, 4> mask_table{1.0, 0.0, 1.0, 1.0};
}  // namespace

namespace snplib {
void CalcAdjustedGRMThread(const uint8_t *geno, size_t num_samples,
                           size_t num_snps, double *covariates,
                           size_t num_covariates, double *matrix,
                           double *gcta_diag) {
  SNP snp(geno, num_samples);
  auto *geno_d = new double[num_samples];
  auto *mask_d = new double[num_samples];
  LogisticRegress<2> worker(num_samples, num_covariates);
  auto m = static_cast<int32_t>(num_samples);
  std::fill(matrix, matrix + num_samples * num_samples, 0.0);
  std::fill(gcta_diag, gcta_diag + num_samples, 0.0);
  for (size_t i = 0; i < num_snps; ++i) {
    snp.UnpackGeno(geno_table, mask_table, geno_d, mask_d);
    worker.Estimate(covariates, geno_d, mask_d);
    auto *u = worker.GetU();
    auto *w = worker.GetW();
    for (size_t j = 0; j < num_samples; ++j) {
      geno_d[j] -= u[j];
      geno_d[j] /= w[j];
      geno_d[j] *= mask_d[j];
    }
    cblas_dsyr(CblasColMajor, CblasLower, m, 1.0, geno_d, 1, matrix, m);
    for (size_t j = 0; j < num_samples; ++j) {
      gcta_diag[j] += (u[j] - 1.0) * geno_d[j] / w[j];
    }
    snp += 1;
  }
  delete[] geno_d;
  delete[] mask_d;
}

void CalcAdjustedGRM(const uint8_t *geno, size_t num_samples, size_t num_snps,
                     double *covariates, size_t num_covariates, double *matrix,
                     double *gcta_diag, size_t num_threads) {
  std::vector<std::thread> workers;
  auto *matrices = new double[num_samples * num_samples * num_threads];
  auto *diags = new double[num_samples * num_threads];
  auto num_snps_job = num_snps / num_threads + 1;
  auto num_snps_left = num_snps % num_threads;
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + (num_samples_left > 0 ? 1 : 0);
  for (size_t i = 0; i < num_snps_left; ++i) {
    workers.emplace_back(CalcAdjustedGRMThread, geno, num_samples, num_snps_job,
                         covariates, num_covariates,
                         matrices + i * num_samples * num_samples,
                         diags + i * num_samples);
    geno += num_snps_job * num_bytes;
  }
  --num_snps_job;
  for (size_t i = num_snps_left; i < num_threads; ++i) {
    workers.emplace_back(CalcAdjustedGRMThread, geno, num_samples, num_snps_job,
                         covariates, num_covariates,
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

void CalcAdmixedGRMThread(const uint8_t *geno, size_t num_samples,
                          size_t num_snps, double *pop_af, double *pop,
                          size_t num_pops, double *matrix, double *gcta_diag) {
  SNP snp(geno, num_samples);
  auto *geno_d = new double[num_samples];
  auto *mask_d = new double[num_samples];
  auto *u = new double[num_samples];
  auto *w = new double[num_samples];
  auto d = static_cast<int32_t>(num_pops);
  auto m = static_cast<int32_t>(num_samples);
  std::fill(matrix, matrix + num_samples * num_samples, 0.0);
  std::fill(gcta_diag, gcta_diag + num_samples, 0.0);
  for (size_t i = 0; i < num_snps; ++i) {
    snp.UnpackGeno(geno_table, mask_table, geno_d, mask_d);
    auto *af = pop_af + i * num_pops;
    cblas_dgemv(CblasColMajor, CblasNoTrans, m, d, 2.0, pop, m, af, 1, 0.0, u,
                1);
    for (size_t j = 0; j < num_samples; ++j) {
      w[j] = std::sqrt(u[j] * (1.0 - u[j] / 2.0));
      geno_d[j] -= u[j];
      geno_d[j] /= w[j];
      geno_d[j] *= mask_d[j];
    }
    cblas_dsyr(CblasColMajor, CblasLower, m, 1.0, geno_d, 1, matrix, m);
    for (size_t j = 0; j < num_samples; ++j) {
      gcta_diag[j] += (u[j] - 1.0) * geno_d[j] / w[j];
    }
  }
  delete[] geno_d;
  delete[] mask_d;
  delete[] u;
  delete[] w;
}

void CalcAdmixedGRM(const uint8_t *geno, size_t num_samples, size_t num_snps,
                    double *pop_af, double *pop, size_t num_pops,
                    double *matrix, double *gcta_diag, size_t num_threads) {
  std::vector<std::thread> workers;
  auto *matrices = new double[num_samples * num_samples * num_threads];
  auto *diags = new double[num_samples * num_threads];
  auto num_snps_job = num_snps / num_threads + 1;
  auto num_snps_left = num_snps % num_threads;
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + (num_samples_left > 0 ? 1 : 0);
  for (size_t i = 0; i < num_snps_left; ++i) {
    workers.emplace_back(CalcAdmixedGRMThread, geno, num_samples, num_snps_job,
                         pop_af, pop, num_pops,
                         matrices + i * num_samples * num_samples,
                         diags + i * num_samples);
    geno += num_snps_job * num_bytes;
    pop_af += num_snps_job * num_pops;
  }
  --num_snps_job;
  for (size_t i = num_snps_left; i < num_threads; ++i) {
    workers.emplace_back(CalcAdmixedGRMThread, geno, num_samples, num_snps_job,
                         pop_af, pop, num_pops,
                         matrices + i * num_samples * num_samples,
                         diags + i * num_samples);
    geno += num_snps_job * num_bytes;
    pop_af += num_snps_job * num_pops;
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
}  // namespace snplib