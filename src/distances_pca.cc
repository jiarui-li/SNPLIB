#include "distances_pca.h"

namespace {
void CalcPairwiseDistances(const double *traits, size_t num_samples,
                           size_t num_traits, size_t num_dims, size_t start_i,
                           size_t start_j, size_t num_pairs,
                           double *distances) {
  auto ii = start_i;
  auto jj = start_j;
  for (size_t l = 0; l < num_pairs; ++l) {
    auto *traits_i = traits + ii * num_dims * num_samples;
    auto *traits_j = traits + jj * num_dims * num_samples;
    auto *dist_ij = distances + l * num_samples;
    for (size_t i = 0; i < num_samples; ++i) {
      double dist = 0.0;
      for (size_t j = 0; j < num_dims; ++j) {
        auto diff =
            traits_i[j * num_samples + i] - traits_j[j * num_samples + i];
        dist += diff * diff;
      }
      dist_ij[i] = std::sqrt(dist);
    }
    jj++;
    if (jj >= num_traits) {
      ii++;
      jj = ii + 1;
    }
  }
}
void CalcDistancesCOVThread() {}
}  // namespace

namespace snplib {
void CalcDistancesPCA(const double *traits, size_t num_samples,
                      size_t num_traits, size_t num_dims, size_t num_components,
                      double *scores, double *vars, double *loadings,
                      size_t num_threads) {}

}  // namespace snplib