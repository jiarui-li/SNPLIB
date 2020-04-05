#ifndef _SNPLIB_SRC_BASIC_STATISTICS_H_
#define _SNPLIB_SRC_BASIC_STATISTICS_H_

#ifdef _MSC_VER
#include <intrin.h>
#else
#include <x86intrin.h>
#endif
#include "snp.h"

namespace snplib {
void CalcAlleleFrequencies(uint8_t *geno, size_t num_samples, size_t num_snps,
                           double *af);
void CalcMissing(uint8_t *geno, size_t num_samples, size_t num_snps,
                 double *ms);
}  // namespace snplib

#endif  //_SNPLIB_SRC_BASIC_STATISTICS_H_
