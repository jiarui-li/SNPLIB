#include "distances_pca.h"

namespace {
std::atomic_size_t ind;
void CalcDistances(const double *traits, size_t num_samples, size_t num_dims,
                   size_t ind_i, size_t ind_j, double *distances) {
  auto *traits_i = traits + ind_i * num_dims * num_samples;
  auto *traits_j = traits + ind_j * num_dims * num_samples;
  std::fill(distances, distances + num_samples, 0.0);
  for (size_t i = 0; i < num_dims; ++i) {
    auto *tmp_i = traits_i + i * num_samples;
    auto *tmp_j = traits_j + i * num_samples;
    for (size_t j = 0; j < num_samples; ++j) {
      distances[j] += (tmp_i[j] - tmp_j[j]) * (tmp_i[j] - tmp_j[j]);
    }
  }
  for (size_t i = 0; i < num_samples; ++i) {
    distances[i] = std::sqrt(distances[i]);
  }
}
void CalcDistancesCOVThread(const double *traits, size_t num_samples,
                            size_t num_traits, size_t num_dims, double *cov) {
  size_t local_ind = ind++;
  size_t local_i = 0;
  size_t local_j = 0;
  while (local_ind > 0) {
    local_ind -= num_traits - local_i;
    local_i++;
  }
  local_i--;
  local_j = num_traits + local_ind - 1;
  std::fill(cov, cov + num_samples * num_samples, 0.0);
  auto *dist = new double[num_samples];
  auto m = static_cast<int32_t>(num_samples);
  while (local_i < num_traits) {
    double mu = std::accumulate(dist, dist + num_samples, 0.0);
    mu /= num_samples;
    for (size_t i = 0; i < num_samples; ++i) {
      dist[i] -= mu;
    }
    cblas_dsyr(CblasColMajor, CblasLower, m, 1.0, dist, 1, cov, m);
    local_ind = ind++;
    local_i = 0;
    while (local_ind > 0) {
      local_ind -= num_traits - local_i;
      local_i++;
    }
    local_i--;
    local_j = num_traits + local_ind - 1;
  }
  delete[] dist;
}
}  // namespace

namespace snplib {
void CalcDistancesPCA(const double *traits, size_t num_samples,
                      size_t num_traits, size_t num_dims, size_t num_components,
                      double *scores, double *vars, double *loadings,
                      size_t num_threads) {
  ind = 1;
}

}  // namespace snplib