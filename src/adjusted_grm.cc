#include "adjusted_grm.h"

namespace {
void CalcAdjustedGRMThread(uint8_t *geno, size_t num_samples, size_t num_snps,
                           double *covariates, size_t num_covariates,
                           double *matrix, double *gcta_diag) {
  snplib::SNP snp(geno, num_samples);
  auto *geno_d = new double[num_samples];
  auto *mask_d = new double[num_samples];
  snplib::LogisticRegress<2> worker(num_samples, num_covariates);
  auto m = static_cast<int32_t>(num_samples);
  std::fill(matrix, matrix + num_samples * num_samples, 0.0);
  std::fill(gcta_diag, gcta_diag + num_samples, 0.0);
  for (size_t l = 0; l < num_snps; ++l) {
    snp.UnpackGeno(geno_d, mask_d);
    worker.Estimate(covariates, geno_d, mask_d);
    auto *u = worker.GetU();
    auto *w = worker.GetW();
    for (size_t i = 0; i < num_samples; ++i) {
      geno_d[i] -= u[i];
      geno_d[i] /= w[i];
      geno_d[i] *= mask_d[i];
    }
    cblas_dsyr(CblasColMajor, CblasLower, m, 1.0, geno_d, 1, matrix, m);
    for (size_t i = 0; i < num_samples; ++i) {
      gcta_diag[i] += (u[i] - 1.0) * geno_d[i] / w[i];
    }
    ++snp;
  }
  delete[] geno_d;
  delete[] mask_d;
}
void CalcAdmixedGRMthread(uint8_t *geno, size_t num_samples, size_t num_snps,
                          double *pop_af, double *pop, size_t num_pops,
                          double *matrix, double *gcta_diag) {
  snplib::SNP snp(geno, num_samples);
  auto *geno_d = new double[num_samples];
  auto *mask_d = new double[num_samples];
  auto *u = new double[num_samples];
  auto *w = new double[num_samples];
  auto d = static_cast<int32_t>(num_pops);
  auto m = static_cast<int32_t>(num_samples);
  std::fill(matrix, matrix + num_samples * num_samples, 0.0);
  std::fill(gcta_diag, gcta_diag + num_samples, 0.0);
  for (size_t l = 0; l < num_snps; ++l) {
    snp.UnpackGeno(geno_d, mask_d);
    auto *af = pop_af + l * num_pops;
    cblas_dgemv(CblasColMajor, CblasNoTrans, m, d, 2.0, pop, m, af, 1, 0.0, u,
                1);
    for (size_t i = 0; i < num_samples; ++i) {
      w[i] = std::sqrt(u[i] * (1.0 - u[i] / 2.0));
      geno_d[i] -= u[i];
      geno_d[i] /= w[i];
      geno_d[i] *= mask_d[i];
    }
    cblas_dsyr(CblasColMajor, CblasLower, m, 1.0, geno_d, 1, matrix, m);
    for (size_t i = 0; i < num_samples; ++i) {
      gcta_diag[i] += (u[i] - 1.0) * geno_d[i] / w[i];
    }
    ++snp;
  }
  delete[] geno_d;
  delete[] mask_d;
  delete[] u;
  delete[] w;
}
}  // namespace

void CalcAdjustedGRM(uint8_t *geno, size_t num_samples, size_t num_snps,
                     double *covariates, size_t num_covariates, double *matrix,
                     double *gcta_diag, size_t num_threads) {
  std::vector<std::thread> workers(num_threads);
  auto *matrices = new double[num_samples * num_samples * num_threads];
  auto *diags = new double[num_samples * num_threads];
  auto num_snps_job = num_snps / num_threads + 1;
  auto num_snps_left = num_snps % num_threads;
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + (num_samples_left > 0 ? 1 : 0);
  for (size_t i = 0; i < num_snps_left; ++i) {
    workers[i] = std::thread(CalcAdjustedGRMThread, geno, num_samples,
                             num_snps_job, covariates, num_covariates,
                             matrices + i * num_samples * num_samples,
                             diags + i * num_samples);
    geno += num_snps_job * num_bytes;
  }
  --num_snps_job;
  for (size_t i = num_snps_left; i < num_threads; ++i) {
    workers[i] = std::thread(CalcAdjustedGRMThread, geno, num_samples,
                             num_snps_job, covariates, num_covariates,
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

void CalcAdmixedGRM(uint8_t *geno, size_t num_samples, size_t num_snps,
                    double *pop_af, double *pop, size_t num_pops,
                    double *matrix, double *gcta_diag, size_t num_threads) {
  std::vector<std::thread> workers(num_threads);
  auto *matrices = new double[num_samples * num_samples * num_threads];
  auto *diags = new double[num_samples * num_threads];
  auto num_snps_job = num_snps / num_threads + 1;
  auto num_snps_left = num_snps % num_threads;
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + (num_samples_left > 0 ? 1 : 0);
  for (size_t i = 0; i < num_snps_left; ++i) {
    workers[i] = std::thread(CalcAdmixedGRMthread, geno, num_samples,
                             num_snps_job, pop_af, pop, num_pops,
                             matrices + i * num_samples * num_samples,
                             diags + i * num_samples);
    pop_af += num_snps_job * num_pops;
    geno += num_snps_job * num_bytes;
  }
  --num_snps_job;
  for (size_t i = num_snps_left; i < num_threads; ++i) {
    workers[i] = std::thread(CalcAdmixedGRMthread, geno, num_samples,
                             num_snps_job, pop_af, pop, num_pops,
                             matrices + i * num_samples * num_samples,
                             diags + i * num_samples);
    pop_af += num_snps_job * num_pops;
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