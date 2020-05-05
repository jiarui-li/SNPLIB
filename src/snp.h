#ifndef _SNPLIB_SRC_SNP_H_
#define _SNPLIB_SRC_SNP_H_

#include <array>
#include <cstring>

namespace snplib {
class SNP {
 private:
  const uint8_t *geno_;
  const size_t num_samples_;
  const size_t num_full_bytes_;
  const size_t num_samples_left_;
  const size_t num_bytes_;

  void ConvertGeno(size_t num_snps, size_t idx, std::array<uint64_t, 4> &g64);
  void ConvertGeno(size_t num_snps, size_t idx, std::array<uint64_t, 5> &g64);

 public:
  SNP(const uint8_t *geno, size_t num_samples);
  SNP(const SNP &) = default;
  SNP(SNP &&) = default;
  SNP &operator=(const SNP &) = delete;
  SNP &operator=(SNP &&) = delete;
  ~SNP() = default;
  uint8_t operator[](size_t idx) const { return geno_[idx]; }
  uint8_t operator()(size_t idx) const {
    auto i = idx / 4;
    auto s = idx % 4;
    return ((*this)[i] >> (2 * s)) & 3u;
  }
  SNP &operator+=(size_t idx) {
    geno_ += idx * num_bytes_;
    return *this;
  }
  template <class T>
  void Copy(T *dest) const {
    memcpy((void *)dest, (const void *)geno_, sizeof(uint8_t) * num_bytes_);
  }
  template <class T>
  void UnpackGeno(const std::array<T, 4> &geno_table, T *geno) {
    for (size_t i = 0; i < num_full_bytes_; ++i) {
      auto t = geno_[i];
      geno[4 * i] = geno_table[t & 3u];
      t >>= 2;
      geno[4 * i + 1] = geno_table[t & 3u];
      t >>= 2;
      geno[4 * i + 2] = geno_table[t & 3u];
      t >>= 2;
      geno[4 * i + 3] = geno_table[t & 3u];
    }
    if (num_samples_left_ > 0u) {
      auto t = geno_[num_full_bytes_];
      for (size_t i = 0; i < num_samples_left_; ++i) {
        geno[4 * num_full_bytes_ + i] = geno_table[t & 3u];
        t >>= 2;
      }
    }
  }
  template <class T>
  void UnpackGeno(const std::array<T, 4> &geno_table,
                  const std::array<T, 4> &mask_table, T *geno, T *mask) {
    for (size_t i = 0; i < num_full_bytes_; ++i) {
      auto t = geno_[i];
      geno[4 * i] = geno_table[t & 3u];
      mask[4 * i] = mask_table[t & 3u];
      t >>= 2;
      geno[4 * i + 1] = geno_table[t & 3u];
      mask[4 * i + 1] = mask_table[t & 3u];
      t >>= 2;
      geno[4 * i + 2] = geno_table[t & 3u];
      mask[4 * i + 2] = mask_table[t & 3u];
      t >>= 2;
      geno[4 * i + 3] = geno_table[t & 3u];
      mask[4 * i + 3] = mask_table[t & 3u];
    }
    if (num_samples_left_ > 0u) {
      auto t = geno_[num_full_bytes_];
      for (size_t i = 0; i < num_samples_left_; ++i) {
        geno[4 * num_full_bytes_ + i] = geno_table[t & 3u];
        mask[4 * num_full_bytes_ + i] = mask_table[t & 3u];
        t >>= 2;
      }
    }
  }
  void TransposeGeno(size_t num_snps, size_t idx, uint64_t *geno64);
  void TransposeGeno(size_t num_snps, uint64_t *geno64);
};
}  // namespace snplib

#endif  //_SNPLIB_SRC_SNP_H_
