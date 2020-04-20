#ifndef _SNPLIB_SRC_SIMULATIONS_H_
#define _SNPLIB_SRC_SIMULATIONS_H_

#include <algorithm>
#include <random>

#include "snp.h"

namespace snplib {
void UpdateAf(const double *aaf, size_t num_pops, size_t num_snps,
              size_t num_generations, size_t effective_sample_size, double *af);
void GenerateIndividuals(const double *af, size_t num_samples, size_t num_snps,
                         uint8_t *geno);
void GenerateAdmixedIndividuals(const double *af, size_t num_snps,
                                size_t num_samples, uint8_t *geno);
void GeneratePairwiseSiblings(uint8_t *parent_geno, size_t num_families,
                              size_t num_snps, uint8_t *siblings_geno);
}  // namespace snplib

#endif  //_SNPLIB_SRC_SIMULATIONS_H_
