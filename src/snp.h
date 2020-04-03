#ifndef SNPLIB_SNP_H
#define SNPLIB_SNP_H

#include <array>
#include <cassert>
#include <cstring>

namespace snplib {

class SNP {
 private:
  uint8_t *geno_;
  const size_t num_samples_;
  const size_t num_full_bytes_;
  const size_t num_samples_left_;
  const size_t num_bytes_;

  void ConvertUGeno(size_t num_snps, size_t idx, std::array<uint64_t, 4> &g64);
  void ConvertSGeno(size_t num_snps, size_t idx, std::array<uint64_t, 5> &g64);

 public:
  SNP(uint8_t *geno, size_t num_samples)
      : geno_(geno),
        num_samples_(num_samples),
        num_full_bytes_(num_samples / 4),
        num_samples_left_(num_samples_ % 4),
        num_bytes_(num_full_bytes_ + (num_samples_left_ != 0 ? 1 : 0)) {}
  SNP(const SNP &) = default;
  SNP(SNP &&) = default;
  SNP &operator=(const SNP &) = delete;
  SNP &operator=(SNP &&) = delete;
  ~SNP() = default;
  uint8_t operator[](size_t idx) const {
    assert(idx < num_samples_);
    auto i = idx / 4;
    auto s = idx % 4;
    return ((*this)[i] >> (2 * s)) & 3u;
  }
  SNP &operator++() {
    geno_ += num_bytes_;
    return *this;
  }
  const SNP operator++(int) {
    SNP result = *this;
    geno_ += num_bytes_;
    return result;
  }
  SNP &operator+=(size_t idx) {
    geno_ += idx * num_bytes_;
    return *this;
  }
  void Unpack(double *geno_d, double *mask_d) const;
  void Unpack(double *geno_d, const std::array<double, 4> &geno_table) const;
  void Flip();
  template <class T>
  void Copy(T *dest) const {
    memcpy((void *)dest, (const void *)geno_, sizeof(uint8_t) * num_bytes_);
  }
  void TransposeUGeno(size_t num_snps, size_t idx, uint64_t *geno64);
  void TransposeSGeno(size_t num_snps, uint64_t *geno64);
};
}  // namespace snplib

#endif  // SNPLIB_SNP_H