#ifndef SNPLIB_BASE_H
#define SNPLIB_BASE_H

#include "snp.h"
#ifdef _MSC_VER
#include <intrin.h>
#else
#include <x86intrin.h>
#endif

extern "C" {
void CalcAlleleFrequencies(uint8_t *geno, size_t num_samples, size_t num_snps,
                           double *af);
void CalcCallrates(uint8_t *geno, size_t num_samples, size_t num_snps,
                   double *cr);
void FlipGeno(uint8_t *geno, size_t num_samples, size_t num_snps,
              const int32_t *idx);
void Keep(uint8_t *src_geno, uint8_t *dest_geno, size_t num_src_samples,
          size_t num_dest_samples, size_t num_snps, const int32_t *idx);
};

#endif  // SNPLIB_BASE_H
