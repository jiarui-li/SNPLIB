#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "data_manage.h"
#include "genetic_variances.h"
#include "gwas.h"
#include "relationships.h"
#include "statistics.h"

namespace py = pybind11;
using array = py::array_t<double, py::array::f_style | py::array::forcecast>;
using geno_t = py::array_t<uint8_t, py::array::f_style | py::array::forcecast>;

// Data Manage
void FlipGeno(geno_t genotype, size_t num_samples, std::vector<int32_t> &idx) {
  auto geno_buf = genotype.request();
  auto *geno_ptr = reinterpret_cast<uint8_t *>(geno_buf.ptr);
  snplib::FlipGeno(geno_ptr, num_samples, idx);
}
void Keep(geno_t src_geno, geno_t dest_geno, size_t num_src_samples,
          size_t num_dest_samples, std::vector<int32_t> &idx) {
  auto src_geno_buf = src_geno.request();
  auto dest_geno_buf = dest_geno.request();
  if (src_geno_buf.shape[1] != dest_geno_buf.shape[1]) {
    throw std::runtime_error("They must have same number of SNPs");
  }
  auto num_snps = src_geno_buf.shape[1];
  auto *src_geno_ptr = reinterpret_cast<uint8_t *>(src_geno_buf.ptr);
  auto *dest_geno_ptr = reinterpret_cast<uint8_t *>(dest_geno_buf.ptr);
  snplib::Keep(src_geno_ptr, dest_geno_ptr, num_src_samples, num_dest_samples,
               num_snps, idx);
}
array UnpackGRMGeno(geno_t genotype, array af, size_t num_samples) {
  auto geno_buf = genotype.request();
  auto *geno_ptr = reinterpret_cast<uint8_t *>(geno_buf.ptr);
  auto num_snps = static_cast<size_t>(geno_buf.shape[1]);
  array geno_d(std::array<size_t, 2>{num_samples, num_snps});
  auto geno_d_buf = geno_d.request();
  auto *geno_d_ptr = reinterpret_cast<double *>(geno_d_buf.ptr);
  auto af_buf = af.request();
  auto *af_ptr = reinterpret_cast<double *>(af_buf.ptr);
  snplib::UnpackGRMGeno(geno_ptr, af_ptr, num_samples, num_snps, geno_d_ptr);
  return geno_d;
}
array UnpackUGeno(geno_t genotype, size_t num_samples) {
  auto geno_buf = genotype.request();
  auto *geno_ptr = reinterpret_cast<uint8_t *>(geno_buf.ptr);
  auto num_snps = static_cast<size_t>(geno_buf.shape[1]);
  array geno_d(std::array<size_t, 2>{num_samples, num_snps});
  auto geno_d_buf = geno_d.request();
  auto *geno_d_ptr = reinterpret_cast<double *>(geno_d_buf.ptr);
  snplib::UnpackUGeno(geno_ptr, num_samples, num_snps, geno_d_ptr);
  return geno_d;
}

// Statistics
array CalcAlleleFrequencies(geno_t genotype, size_t num_samples) {
  auto geno_buf = genotype.request();
  auto *geno_ptr = reinterpret_cast<uint8_t *>(geno_buf.ptr);
  auto num_snps = static_cast<size_t>(geno_buf.shape[1]);
  array af(num_snps);
  auto af_buf = af.request();
  auto *af_ptr = reinterpret_cast<double *>(af_buf.ptr);
  snplib::CalcAlleleFrequencies(geno_ptr, num_samples, num_snps, af_ptr);
  return af;
}
array CalcMissing(geno_t genotype, size_t num_samples) {
  auto geno_buf = genotype.request();
  auto *geno_ptr = reinterpret_cast<uint8_t *>(geno_buf.ptr);
  auto num_snps = static_cast<size_t>(geno_buf.shape[1]);
  array ms(num_snps);
  auto ms_buf = ms.request();
  auto *ms_ptr = reinterpret_cast<double *>(ms_buf.ptr);
  snplib::CalcMissing(geno_ptr, num_samples, num_snps, ms_ptr);
  return ms;
}
array CalcAdjustedAF(geno_t genotype, array covariates, size_t num_threads) {
  auto geno_buf = genotype.request();
  auto *geno_ptr = reinterpret_cast<uint8_t *>(geno_buf.ptr);
  auto num_snps = static_cast<size_t>(geno_buf.shape[1]);
  auto cov_buf = covariates.request();
  size_t num_samples = cov_buf.shape[0];
  auto num_covariates = cov_buf.shape[1];
  auto *cov_ptr = reinterpret_cast<double *>(cov_buf.ptr);
  array af(std::array<size_t, 2>{num_samples, num_snps});
  auto af_buf = af.request();
  auto *af_ptr = reinterpret_cast<double *>(af_buf.ptr);
  snplib::CalcAdjustedAF(geno_ptr, num_samples, num_snps, cov_ptr,
                         num_covariates, af_ptr, num_threads);
  return af;
}
array CalcAdjustedMAF(geno_t genotype, array covariates, size_t num_threads) {
  auto geno_buf = genotype.request();
  auto *geno_ptr = reinterpret_cast<uint8_t *>(geno_buf.ptr);
  auto num_snps = static_cast<size_t>(geno_buf.shape[1]);
  auto cov_buf = covariates.request();
  auto num_samples = cov_buf.shape[0];
  auto num_covariates = cov_buf.shape[1];
  auto *cov_ptr = reinterpret_cast<double *>(cov_buf.ptr);
  array af(num_snps);
  auto af_buf = af.request();
  auto *af_ptr = reinterpret_cast<double *>(af_buf.ptr);
  snplib::CalcAdjustedMAF(geno_ptr, num_samples, num_snps, cov_ptr,
                          num_covariates, af_ptr, num_threads);
  return af;
}

// Relationships
std::tuple<array, array> CalcAdjustedGRM(geno_t genotype, array covariates,
                                         size_t num_threads) {
  auto geno_buf = genotype.request();
  auto *geno_ptr = reinterpret_cast<uint8_t *>(geno_buf.ptr);
  auto num_snps = static_cast<size_t>(geno_buf.shape[1]);
  auto cov_buf = covariates.request();
  size_t num_samples = cov_buf.shape[0];
  auto num_covariates = cov_buf.shape[1];
  auto *cov_ptr = reinterpret_cast<double *>(cov_buf.ptr);
  array matrix(std::array<size_t, 2>{num_samples, num_samples});
  array gcta_diag(num_samples);
  auto matrix_buf = matrix.request();
  auto *matrix_ptr = reinterpret_cast<double *>(matrix_buf.ptr);
  auto gcta_buf = gcta_diag.request();
  auto *gcta_ptr = reinterpret_cast<double *>(gcta_buf.ptr);
  snplib::CalcAdjustedGRM(geno_ptr, num_samples, num_snps, cov_ptr,
                          num_covariates, matrix_ptr, gcta_ptr, num_threads);
  return std::make_tuple(matrix, gcta_diag);
}
std::tuple<array, array> CalcAdmixedGRM(geno_t genotype, array pop_af,
                                        array pop, size_t num_threads) {
  auto geno_buf = genotype.request();
  auto *geno_ptr = reinterpret_cast<uint8_t *>(geno_buf.ptr);
  auto num_snps = static_cast<size_t>(geno_buf.shape[1]);
  auto pop_af_buf = pop_af.request();
  auto num_pops = static_cast<size_t>(pop_af_buf.shape[0]);
  auto *pop_af_ptr = reinterpret_cast<double *>(pop_af_buf.ptr);
  auto pop_buf = pop.request();
  auto *pop_ptr = reinterpret_cast<double *>(pop_buf.ptr);
  auto num_samples = static_cast<size_t>(pop_buf.shape[0]);
  array matrix(std::array<size_t, 2>{num_samples, num_samples});
  array gcta_diag(num_samples);
  auto matrix_buf = matrix.request();
  auto *matrix_ptr = reinterpret_cast<double *>(matrix_buf.ptr);
  auto gcta_buf = gcta_diag.request();
  auto *gcta_ptr = reinterpret_cast<double *>(gcta_buf.ptr);
  snplib::CalcAdmixedGRM(geno_ptr, num_samples, num_snps, pop_af_ptr, pop_ptr,
                         num_pops, matrix_ptr, gcta_ptr, num_threads);
  return std::make_tuple(matrix, gcta_diag);
}
array CalcGRMMatrix(geno_t genotype, array af, size_t num_samples,
                    size_t num_threads) {
  auto geno_buf = genotype.request();
  auto *geno_ptr = reinterpret_cast<uint8_t *>(geno_buf.ptr);
  auto num_snps = static_cast<size_t>(geno_buf.shape[1]);
  array matrix(std::array<size_t, 2>{num_samples, num_samples});
  auto matrix_buf = matrix.request();
  auto *matrix_ptr = reinterpret_cast<double *>(matrix_buf.ptr);
  auto af_buf = af.request();
  auto *af_ptr = reinterpret_cast<double *>(af_buf.ptr);
  snplib::CalcGRMMatrix(geno_ptr, af_ptr, num_samples, num_snps, matrix_ptr,
                        num_threads);
  return matrix;
}
array CalcGCTADiagonal(geno_t genotype, array af, size_t num_samples,
                       size_t num_threads) {
  auto geno_buf = genotype.request();
  auto *geno_ptr = reinterpret_cast<uint8_t *>(geno_buf.ptr);
  auto num_snps = static_cast<size_t>(geno_buf.shape[1]);
  array gcta_diag(num_samples);
  auto gcta_buf = gcta_diag.request();
  auto *gcta_ptr = reinterpret_cast<double *>(gcta_buf.ptr);
  auto af_buf = af.request();
  auto *af_ptr = reinterpret_cast<double *>(af_buf.ptr);
  snplib::CalcGCTADiagonal(geno_ptr, af_ptr, num_samples, num_snps, gcta_ptr,
                           num_threads);
  return gcta_diag;
}
array CalcIBSMatrix(geno_t genotype, size_t num_samples, size_t num_threads) {
  auto geno_buf = genotype.request();
  auto *geno_ptr = reinterpret_cast<uint8_t *>(geno_buf.ptr);
  auto num_snps = static_cast<size_t>(geno_buf.shape[1]);
  array matrix(std::array<size_t, 2>{num_samples, num_samples});
  array gcta_diag(num_samples);
  auto matrix_buf = matrix.request();
  auto *matrix_ptr = reinterpret_cast<double *>(matrix_buf.ptr);
  snplib::CalcIBSMatrix(geno_ptr, num_samples, num_snps, matrix_ptr,
                        num_threads);
  return matrix;
}
array CalcIBSConnection(geno_t src_geno, geno_t dest_geno,
                        size_t num_src_samples, size_t num_dest_samples,
                        size_t num_threads) {
  auto src_geno_buf = src_geno.request();
  auto dest_geno_buf = dest_geno.request();
  if (src_geno_buf.shape[1] != dest_geno_buf.shape[1]) {
    throw std::runtime_error("They must have same number of SNPs");
  }
  auto num_snps = src_geno_buf.shape[1];
  auto *src_geno_ptr = reinterpret_cast<uint8_t *>(src_geno_buf.ptr);
  auto *dest_geno_ptr = reinterpret_cast<uint8_t *>(dest_geno_buf.ptr);
  array connection(num_dest_samples);
  auto connect_buf = connection.request();
  auto *connect_ptr = reinterpret_cast<double *>(connect_buf.ptr);
  snplib::CalcIBSConnection(src_geno_ptr, num_src_samples, dest_geno_ptr,
                            num_dest_samples, num_snps, connect_ptr,
                            num_threads);
  return connection;
}
array CalcKINGMatrix(geno_t genotype, size_t num_samples, size_t num_threads) {
  auto geno_buf = genotype.request();
  auto *geno_ptr = reinterpret_cast<uint8_t *>(geno_buf.ptr);
  auto num_snps = static_cast<size_t>(geno_buf.shape[1]);
  array matrix(std::array<size_t, 2>{num_samples, num_samples});
  array gcta_diag(num_samples);
  auto matrix_buf = matrix.request();
  auto *matrix_ptr = reinterpret_cast<double *>(matrix_buf.ptr);
  snplib::CalcKINGMatrix(geno_ptr, num_samples, num_snps, matrix_ptr,
                         num_threads);
  return matrix;
}
array CalcUGRMMatrix(geno_t genotype, size_t num_samples, size_t num_threads) {
  auto geno_buf = genotype.request();
  auto *geno_ptr = reinterpret_cast<uint8_t *>(geno_buf.ptr);
  auto num_snps = static_cast<size_t>(geno_buf.shape[1]);
  array matrix(std::array<size_t, 2>{num_samples, num_samples});
  array gcta_diag(num_samples);
  auto matrix_buf = matrix.request();
  auto *matrix_ptr = reinterpret_cast<double *>(matrix_buf.ptr);
  snplib::CalcUGRMMatrix(geno_ptr, num_samples, num_snps, matrix_ptr,
                         num_threads);
  return matrix;
}

// Genetic Variances
std::tuple<array, array> CalcUniLMM(array traits, array covariates,
                                    array lambda, size_t num_threads) {
  auto traits_buf = traits.request();
  auto num_samples = static_cast<size_t>(traits_buf.shape[0]);
  auto num_traits = static_cast<size_t>(traits_buf.shape[1]);
  auto *traits_ptr = reinterpret_cast<double *>(traits_buf.ptr);
  auto cov_buf = covariates.request();
  auto num_covariates = static_cast<size_t>(cov_buf.shape[1]);
  auto *cov_ptr = reinterpret_cast<double *>(cov_buf.ptr);
  auto lambda_buf = lambda.request();
  auto *lambda_ptr = reinterpret_cast<double *>(lambda_buf.ptr);
  array vars(std::array<size_t, 2>{2, num_traits});
  array res(std::array<size_t, 2>{num_samples, num_traits});
  auto vars_buf = vars.request();
  auto *vars_ptr = reinterpret_cast<double *>(vars_buf.ptr);
  auto res_buf = res.request();
  auto *res_ptr = reinterpret_cast<double *>(res_buf.ptr);
  snplib::CalcUniLMM(traits_ptr, cov_ptr, lambda_ptr, num_samples,
                     num_covariates, num_traits, vars_ptr, res_ptr,
                     num_threads);
  return std::make_tuple(vars, res);
}
std::tuple<array, array> CalcMLMMSigmas(array traits, array covariates,
                                        array lambda, array res, array vars,
                                        size_t num_dims, size_t num_threads) {
  auto traits_buf = traits.request();
  auto num_samples = static_cast<size_t>(traits_buf.shape[0]);
  auto num_traits = static_cast<size_t>(traits_buf.shape[1]) / num_dims;
  auto *traits_ptr = reinterpret_cast<double *>(traits_buf.ptr);
  auto cov_buf = covariates.request();
  auto num_covariates = static_cast<size_t>(cov_buf.shape[1]);
  auto *cov_ptr = reinterpret_cast<double *>(cov_buf.ptr);
  auto lambda_buf = lambda.request();
  auto *lambda_ptr = reinterpret_cast<double *>(lambda_buf.ptr);
  auto res_buf = res.request();
  auto *res_ptr = reinterpret_cast<double *>(res_buf.ptr);
  auto vars_buf = vars.request();
  auto *vars_ptr = reinterpret_cast<double *>(vars_buf.ptr);
  array sigma_e(std::array<size_t, 3>{num_dims, num_dims, num_traits});
  array sigma_g(std::array<size_t, 3>{num_dims, num_dims, num_traits});
  auto sigma_e_buf = sigma_e.request();
  auto *sigma_e_ptr = reinterpret_cast<double *>(sigma_e_buf.ptr);
  auto sigma_g_buf = sigma_g.request();
  auto *sigma_g_ptr = reinterpret_cast<double *>(sigma_g_buf.ptr);
  snplib::CalcMLMMSigmas(cov_ptr, lambda_ptr, traits_ptr, res_ptr, vars_ptr,
                         num_samples, num_covariates, num_dims, num_traits,
                         sigma_e_ptr, sigma_g_ptr, num_threads);
  return std::make_tuple(sigma_e, sigma_g);
}
std::tuple<array, array> CalcRMLMMSigmas(array lambda, array res, array vars,
                                         size_t num_covariates, size_t num_dims,
                                         size_t num_threads) {
  auto lambda_buf = lambda.request();
  auto *lambda_ptr = reinterpret_cast<double *>(lambda_buf.ptr);
  auto res_buf = res.request();
  auto num_samples = static_cast<size_t>(res_buf.shape[0]);
  auto num_traits = static_cast<size_t>(res_buf.shape[1]) / num_dims;
  auto *res_ptr = reinterpret_cast<double *>(res_buf.ptr);
  auto vars_buf = vars.request();
  auto *vars_ptr = reinterpret_cast<double *>(vars_buf.ptr);
  array sigma_e(std::array<size_t, 3>{num_dims, num_dims, num_traits});
  array sigma_g(std::array<size_t, 3>{num_dims, num_dims, num_traits});
  auto sigma_e_buf = sigma_e.request();
  auto *sigma_e_ptr = reinterpret_cast<double *>(sigma_e_buf.ptr);
  auto sigma_g_buf = sigma_g.request();
  auto *sigma_g_ptr = reinterpret_cast<double *>(sigma_g_buf.ptr);
  snplib::CalcRMLMMSigmas(lambda_ptr, res_ptr, vars_ptr, num_samples,
                          num_covariates, num_dims, num_traits, sigma_e_ptr,
                          sigma_g_ptr, num_threads);
  return std::make_tuple(sigma_e, sigma_g);
}

// GWAS TODO