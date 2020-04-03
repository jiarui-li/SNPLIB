#include "logistic_utils.h"

namespace {
void UnpackResGenoThread(uint8_t *geno, size_t num_samples, size_t num_snps,
                         double *covariates, size_t num_covariates,
                         double *geno_d) {
  snplib::SNP snp(geno, num_samples);
  auto *mask_d = new double[num_samples];
  snplib::LogisticRegress<2> worker(num_samples, num_covariates);
  for (size_t l = 0; l < num_snps; ++l) {
    auto *tmp_g = geno_d + l * num_samples;
    snp.UnpackGeno(tmp_g, mask_d);
    worker.Estimate(covariates, tmp_g, mask_d);
    auto *u = worker.GetU();
    auto *w = worker.GetW();
    for (size_t i = 0; i < num_samples; ++i) {
      tmp_g[i] -= u[i];
      tmp_g[i] /= w[i];
      tmp_g[i] *= mask_d[i];
    }
    ++snp;
  }
  delete[] mask_d;
}
void CalcLogisticGWASThread(const double *trait, const double *covariates,
                            uint8_t *geno, size_t num_samples,
                            size_t num_covariates, size_t num_snps,
                            double *chi2stat) {
  auto *cov = new double[num_samples * (num_covariates + 1)];
  auto *geno_d = cov + num_samples * num_covariates;
  auto *mask_d = new double[num_samples];
  std::copy(covariates, covariates + num_samples * num_covariates, cov);
  snplib::SNP snp(geno, num_samples);
  snplib::LogisticRegress<1> worker_full(num_samples, num_covariates + 1);
  snplib::LogisticRegress<1> worker_reduce(num_samples, num_covariates);
  auto m = static_cast<int32_t>(num_samples);
  for (size_t l = 0; l < num_snps; ++l) {
    snp.UnpackGeno(geno_d, mask_d);
    worker_full.Estimate(cov, trait, mask_d);
    worker_reduce.Estimate(covariates, trait, mask_d);
    chi2stat[l] = worker_full.CalcLikelihood(cov, trait, mask_d) -
                  worker_reduce.CalcLikelihood(covariates, trait, mask_d);
    ++snp;
  }
  delete[] mask_d;
}
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
void CalcSNPGFTThread(const double *covariates, const double *af, uint8_t *geno,
                      size_t num_samples, size_t num_covariates,
                      size_t num_snps, double *gft) {
  snplib::SNP snp(geno, num_samples);
  auto *geno_d = new double[num_samples];
  auto *mask_d = new double[num_samples];
  snplib::LogisticRegress<2> worker(num_samples, num_covariates);
  for (size_t l = 0; l < num_snps; ++l) {
    auto aff = std::log(af[l]);
    auto laff = std::log(1.0 - af[l]);
    snp.UnpackGeno(geno_d, mask_d);
    worker.Estimate(covariates, geno_d, mask_d);
    gft[l] = worker.CalcLikelihood(covariates, geno_d, mask_d);
    double gg = 0.0;
    for (size_t i = 0; i < num_samples; ++i) {
      gg += mask_d[i] * (geno_d[i] * aff + (2.0 - geno_d[i]) * laff);
    }
    gft[l] -= gg;
    ++snp;
  }
  delete[] geno_d;
  delete[] mask_d;
}

}  // namespace

void UnpackResGeno(uint8_t *geno, size_t num_samples, size_t num_snps,
                   double *covariates, size_t num_covariates, double *geno_d,
                   size_t num_threads) {
  std::vector<std::thread> workers(num_threads);
  auto num_snps_job = num_snps / num_threads + 1;
  auto num_snps_left = num_snps % num_threads;
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + (num_samples_left > 0 ? 1 : 0);
  for (size_t i = 0; i < num_snps_left; ++i) {
    workers[i] = std::thread(UnpackResGenoThread, geno, num_samples,
                             num_snps_job, covariates, num_covariates, geno_d);
    geno += num_bytes * num_snps_job;
    geno_d += num_samples * num_snps_job;
  }
  --num_snps_job;
  for (size_t i = num_snps_left; i < num_threads; ++i) {
    workers[i] = std::thread(UnpackResGenoThread, geno, num_samples,
                             num_snps_job, covariates, num_covariates, geno_d);
    geno += num_bytes * num_snps_job;
    geno_d += num_samples * num_snps_job;
  }
  for (auto &&iter : workers) {
    iter.join();
  }
}

void CalcLogisticGWAS(const double *trait, const double *covariates,
                      uint8_t *geno, size_t num_samples, size_t num_covariates,
                      size_t num_snps, double *chi2stat, size_t num_threads) {
  std::vector<std::thread> workers(num_threads);
  auto num_snps_job = num_snps / num_threads + 1;
  auto num_snps_left = num_snps % num_threads;
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + (num_samples_left > 0 ? 1 : 0);
  for (size_t i = 0; i < num_snps_left; ++i) {
    workers[i] =
        std::thread(CalcLogisticGWASThread, trait, covariates, geno,
                    num_samples, num_covariates, num_snps_job, chi2stat);
    geno += num_bytes * num_snps_job;
    chi2stat += num_snps_job;
  }
  --num_snps_job;
  for (size_t i = num_snps_left; i < num_threads; ++i) {
    workers[i] =
        std::thread(CalcLogisticGWASThread, trait, covariates, geno,
                    num_samples, num_covariates, num_snps_job, chi2stat);
    geno += num_bytes * num_snps_job;
    chi2stat += num_snps_job;
  }
  for (auto &&iter : workers) {
    iter.join();
  }
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

void CalcSNPGFT(const double *covariates, const double *af, uint8_t *geno,
                size_t num_samples, size_t num_covariates, size_t num_snps,
                double *gft, size_t num_threads) {
  std::vector<std::thread> workers(num_threads);
  auto num_snps_job = num_snps / num_threads + 1;
  auto num_snps_left = num_snps % num_threads;
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + (num_samples_left > 0 ? 1 : 0);
  for (size_t i = 0; i < num_snps_left; ++i) {
    workers[i] = std::thread(CalcSNPGFTThread, covariates, af, geno,
                             num_samples, num_covariates, num_snps_job, gft);
    geno += num_bytes * num_snps_job;
    af += num_snps_job;
    gft += num_snps_job;
  }
  --num_snps_job;
  for (size_t i = num_snps_left; i < num_threads; ++i) {
    workers[i] = std::thread(CalcSNPGFTThread, covariates, af, geno,
                             num_samples, num_covariates, num_snps_job, gft);
    geno += num_bytes * num_snps_job;
    af += num_snps_job;
    gft += num_snps_job;
  }
  for (auto &&iter : workers) {
    iter.join();
  }
}