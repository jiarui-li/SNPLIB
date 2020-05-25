#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "adjusted_grm.h"
#include "data_manage.h"
#include "genetic_variances.h"
#include "grm.h"
#include "gwas.h"
#include "ibs.h"
#include "king.h"
#include "simulations.h"
#include "statistics.h"
#include "ugrm.h"

namespace py = pybind11;
using array = py::array_t<double, py::array::f_style | py::array::forcecast>;
using geno_t = py::array_t<uint8_t, py::array::f_style | py::array::forcecast>;
using index_t = py::array_t<int32_t, py::array::f_style | py::array::forcecast>;

// Data Manage
void FlipGeno(geno_t genotype, size_t num_samples, index_t index) {
  auto geno_buf = genotype.request();
  auto num_snps = static_cast<size_t>(geno_buf.shape[1]);
  auto *geno_ptr = reinterpret_cast<uint8_t *>(geno_buf.ptr);
  auto idx_buf = index.request();
  auto *idx_ptr = reinterpret_cast<int32_t *>(idx_buf.ptr);
  snplib::FlipGeno(geno_ptr, num_samples, num_snps, idx_ptr);
}
geno_t Keep(geno_t geno, size_t num_samples, index_t index) {
  auto geno_buf = geno.request();
  auto num_snps = static_cast<size_t>(geno_buf.shape[1]);
  auto *geno_ptr = reinterpret_cast<uint8_t *>(geno_buf.ptr);
  auto idx_buf = index.request();
  auto *idx_ptr = reinterpret_cast<int32_t *>(idx_buf.ptr);
  auto num_dest_samples = static_cast<size_t>(idx_buf.size);
  auto num_bytes = num_dest_samples / 4 + ((num_dest_samples % 4) > 0 ? 1 : 0);
  auto dest_geno = geno_t(std::array<size_t, 2>{num_bytes, num_snps});
  auto dest_geno_buf = dest_geno.request();
  auto *dest_geno_ptr = reinterpret_cast<uint8_t *>(dest_geno_buf.ptr);
  snplib::Keep(geno_ptr, dest_geno_ptr, num_samples, num_dest_samples, num_snps,
               idx_ptr);
  return dest_geno;
}
array UnpackGeno(geno_t genotype, size_t num_samples) {
  auto geno_buf = genotype.request();
  auto *geno_ptr = reinterpret_cast<uint8_t *>(geno_buf.ptr);
  auto num_snps = static_cast<size_t>(geno_buf.shape[1]);
  array geno_d(std::array<size_t, 2>{num_samples, num_snps});
  auto geno_d_buf = geno_d.request();
  auto *geno_d_ptr = reinterpret_cast<double *>(geno_d_buf.ptr);
  snplib::UnpackGeno(geno_ptr, num_samples, num_snps, geno_d_ptr);
  return geno_d;
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
index_t FindUnrelatedGroup(array matrix, double threshold) {
  auto matrix_buf = matrix.request();
  auto *matrix_ptr = reinterpret_cast<double *>(matrix_buf.ptr);
  auto num_samples = static_cast<size_t>(matrix_buf.shape[0]);
  auto idx = snplib::FindUnrelatedGroup(matrix_ptr, num_samples, threshold);
  index_t result(idx.size());
  auto result_buf = result.request();
  auto *result_ptr = reinterpret_cast<int32_t *>(result_buf.ptr);
  std::copy(idx.begin(), idx.end(), result_ptr);
  return result;
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
std::tuple<array, array> CalcMultiLMM_REML(array traits, array covariates,
                                           array lambda, array res, array vars,
                                           size_t num_dims,
                                           size_t num_threads) {
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
  snplib::CalcMultiLMM_REML(cov_ptr, lambda_ptr, traits_ptr, res_ptr, vars_ptr,
                            num_samples, num_covariates, num_dims, num_traits,
                            sigma_e_ptr, sigma_g_ptr, num_threads);
  return std::make_tuple(sigma_e, sigma_g);
}
std::tuple<array, array> CalcMultiLMM_RML(array lambda, array res, array vars,
                                          size_t num_covariates,
                                          size_t num_dims, size_t num_threads) {
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
  snplib::CalcMultiLMM_RML(lambda_ptr, res_ptr, vars_ptr, num_samples,
                           num_covariates, num_dims, num_traits, sigma_e_ptr,
                           sigma_g_ptr, num_threads);
  return std::make_tuple(sigma_e, sigma_g);
}

// GWAS
std::tuple<array, array> CalcLinearRegressionGWAS(geno_t genotype,
                                                  array covariates, array trait,
                                                  size_t num_threads) {
  auto geno_buf = genotype.request();
  auto *geno_ptr = reinterpret_cast<uint8_t *>(geno_buf.ptr);
  auto num_snps = static_cast<size_t>(geno_buf.shape[1]);
  auto cov_buf = covariates.request();
  auto num_samples = static_cast<size_t>(cov_buf.shape[0]);
  auto num_covariates = static_cast<size_t>(cov_buf.shape[1]);
  auto *cov_ptr = reinterpret_cast<double *>(cov_buf.ptr);
  auto trait_buf = trait.request();
  auto *trait_ptr = reinterpret_cast<double *>(trait_buf.ptr);
  array betas(num_snps);
  array stats(num_snps);
  auto betas_buf = betas.request();
  auto *betas_ptr = reinterpret_cast<double *>(betas_buf.ptr);
  auto stats_buf = stats.request();
  auto *stats_ptr = reinterpret_cast<double *>(stats_buf.ptr);
  snplib::CalcLinearRegressionGWAS(geno_ptr, num_samples, num_snps, cov_ptr,
                                   num_covariates, trait_ptr, betas_ptr,
                                   stats_ptr, num_threads);
  return std::make_tuple(betas, stats);
}

std::tuple<array, array> CalcLogisticGWAS(geno_t genotype, array covariates,
                                          array trait, size_t num_threads) {
  auto geno_buf = genotype.request();
  auto *geno_ptr = reinterpret_cast<uint8_t *>(geno_buf.ptr);
  auto num_snps = static_cast<size_t>(geno_buf.shape[1]);
  auto cov_buf = covariates.request();
  auto num_samples = static_cast<size_t>(cov_buf.shape[0]);
  auto num_covariates = static_cast<size_t>(cov_buf.shape[1]);
  auto *cov_ptr = reinterpret_cast<double *>(cov_buf.ptr);
  auto trait_buf = trait.request();
  auto *trait_ptr = reinterpret_cast<double *>(trait_buf.ptr);
  array betas(num_snps);
  array stats(num_snps);
  auto betas_buf = betas.request();
  auto *betas_ptr = reinterpret_cast<double *>(betas_buf.ptr);
  auto stats_buf = stats.request();
  auto *stats_ptr = reinterpret_cast<double *>(stats_buf.ptr);
  snplib::CalcLogisticGWAS(geno_ptr, num_samples, num_snps, cov_ptr,
                           num_covariates, trait_ptr, betas_ptr, stats_ptr,
                           num_threads);
  return std::make_tuple(betas, stats);
}

std::tuple<array, array> CalcCCAGWAS(geno_t genotype, array trait,
                                     size_t num_threads) {
  auto geno_buf = genotype.request();
  auto *geno_ptr = reinterpret_cast<uint8_t *>(geno_buf.ptr);
  auto num_snps = static_cast<size_t>(geno_buf.shape[1]);
  auto trait_buf = trait.request();
  auto num_samples = static_cast<size_t>(trait_buf.shape[0]);
  auto num_dims = static_cast<size_t>(trait_buf.shape[1]);
  auto *trait_ptr = reinterpret_cast<double *>(trait_buf.ptr);
  array betas(std::array<size_t, 2>{num_dims, num_snps});
  array rho2(num_snps);
  auto betas_buf = betas.request();
  auto *betas_ptr = reinterpret_cast<double *>(betas_buf.ptr);
  auto rho2_buf = rho2.request();
  auto *rho2_ptr = reinterpret_cast<double *>(rho2_buf.ptr);
  snplib::CalcCCAGWAS(geno_ptr, num_samples, num_snps, trait_ptr, num_dims,
                      betas_ptr, rho2_ptr, num_threads);
  return std::make_tuple(betas, rho2);
}
array CalcCCAReplication(geno_t genotype, array scores, array betas,
                         size_t num_threads) {
  auto geno_buf = genotype.request();
  auto *geno_ptr = reinterpret_cast<uint8_t *>(geno_buf.ptr);
  auto num_snps = static_cast<size_t>(geno_buf.shape[1]);
  auto scores_buf = scores.request();
}
std::tuple<array, array, array> CalcUniLMMGWAS(geno_t genotype,
                                               array covariates, array trait,
                                               array lambda, array V,
                                               size_t num_threads) {
  auto geno_buf = genotype.request();
  auto *geno_ptr = reinterpret_cast<uint8_t *>(geno_buf.ptr);
  auto num_snps = static_cast<size_t>(geno_buf.shape[1]);
  auto cov_buf = covariates.request();
  auto num_samples = static_cast<size_t>(cov_buf.shape[0]);
  auto num_covariates = static_cast<size_t>(cov_buf.shape[1]);
  auto *cov_ptr = reinterpret_cast<double *>(cov_buf.ptr);
  auto trait_buf = trait.request();
  auto *trait_ptr = reinterpret_cast<double *>(trait_buf.ptr);
  auto lambda_buf = lambda.request();
  auto *lambda_ptr = reinterpret_cast<double *>(lambda_buf.ptr);
  auto V_buf = V.request();
  auto *V_ptr = reinterpret_cast<double *>(V_buf.ptr);
  array betas(num_snps);
  array fstats(num_snps);
  array dfs(num_snps);
  auto betas_buf = betas.request();
  auto *betas_ptr = reinterpret_cast<double *>(betas_buf.ptr);
  auto fstats_buf = fstats.request();
  auto *fstats_ptr = reinterpret_cast<double *>(fstats_buf.ptr);
  auto dfs_buf = dfs.request();
  auto *dfs_ptr = reinterpret_cast<double *>(dfs_buf.ptr);
  snplib::CalcUniLMMGWAS(geno_ptr, num_samples, num_snps, lambda_ptr, V_ptr,
                         cov_ptr, num_covariates, trait_ptr, betas_ptr,
                         fstats_ptr, dfs_ptr, num_threads);
  return std::make_tuple(betas, fstats, dfs);
}
// Simulations
array UpdateAf(array aaf, size_t num_pops, size_t num_generations,
               size_t effective_sample_size) {
  auto aaf_buf = aaf.request();
  auto *aaf_ptr = reinterpret_cast<double *>(aaf_buf.ptr);
  auto num_snps = static_cast<size_t>(aaf_buf.size);
  array af(std::array<size_t, 2>{num_snps, num_pops});
  auto af_buf = af.request();
  auto *af_ptr = reinterpret_cast<double *>(af_buf.ptr);
  snplib::UpdateAf(aaf_ptr, num_pops, num_snps, num_generations,
                   effective_sample_size, af_ptr);
  return af;
}
geno_t GenerateIndividuals(array af) {
  auto af_buf = af.request();
  auto *af_ptr = reinterpret_cast<double *>(af_buf.ptr);
  auto num_samples = static_cast<size_t>(af_buf.shape[0]);
  auto num_snps = static_cast<size_t>(af_buf.shape[1]);
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + (num_samples_left != 0 ? 1 : 0);
  geno_t geno(std::array<size_t, 2>{num_bytes, num_snps});
  auto geno_buf = geno.request();
  auto *geno_ptr = reinterpret_cast<uint8_t *>(geno_buf.ptr);
  snplib::GenerateIndividuals(af_ptr, num_samples, num_snps, geno_ptr);
  return geno;
}
geno_t GenerateAdmixedIndividuals(array af, size_t num_samples) {
  auto af_buf = af.request();
  auto *af_ptr = reinterpret_cast<double *>(af_buf.ptr);
  auto num_snps = static_cast<size_t>(af_buf.shape[1]);
  auto num_full_bytes = num_samples / 4;
  auto num_samples_left = num_samples % 4;
  auto num_bytes = num_full_bytes + (num_samples_left != 0 ? 1 : 0);
  geno_t geno(std::array<size_t, 2>{num_bytes, num_snps});
  auto geno_buf = geno.request();
  auto *geno_ptr = reinterpret_cast<uint8_t *>(geno_buf.ptr);
  snplib::GenerateAdmixedIndividuals(af_ptr, num_snps, num_samples, geno_ptr);
  return geno;
}
geno_t GeneratePairwiseSiblings(geno_t parent_geno, size_t num_families) {
  auto geno_p_buf = parent_geno.request();
  auto num_snps = static_cast<size_t>(geno_p_buf.shape[1]);
  auto *geno_p_ptr = reinterpret_cast<uint8_t *>(geno_p_buf.ptr);
  auto num_bytes = static_cast<size_t>(geno_p_buf.shape[0]);
  geno_t geno(std::array<size_t, 2>{num_bytes, num_snps});
  auto geno_buf = geno.request();
  auto *geno_ptr = reinterpret_cast<uint8_t *>(geno_buf.ptr);
  snplib::GeneratePairwiseSiblings(geno_p_ptr, num_families, num_snps,
                                   geno_ptr);
  return geno;
}
// Pybind11 interface
PYBIND11_MODULE(_SNPLIB, m) {
  m.def("FlipGeno", &FlipGeno, "Flip the genotype");
  m.def("Keep", &Keep, "Keep the selected individuals");
  m.def("UnpackGeno", &UnpackGeno, "");
  m.def("UnpackGRMGeno", &UnpackGRMGeno,
        "Unpack the genotype data into the standardized format");
  m.def("UnpackUGeno", &UnpackUGeno,
        "Unpack the genotype data into the unnormalized format");
  m.def("CalcAlleleFrequencies", &CalcAlleleFrequencies,
        "Calculate SNPs' allele frequencies");
  m.def("CalcMissing", &CalcMissing, "Calculate SNPs' missingness");
  m.def("CalcAdjustedAF", &CalcAdjustedAF,
        "Calculate SNPs' ancestry adjusted AF");
  m.def("CalcAdjustedMAF", &CalcAdjustedMAF, "");
  m.def("CalcAdjustedGRM", &CalcAdjustedGRM, "");
  m.def("CalcAdmixedGRM", &CalcAdmixedGRM, "");
  m.def("CalcGRMMatrix", &CalcGRMMatrix, "");
  m.def("CalcGCTADiagonal", &CalcGCTADiagonal, "");
  m.def("CalcIBSMatrix", &CalcIBSMatrix, "");
  m.def("CalcIBSConnection", &CalcIBSConnection, "");
  m.def("CalcKINGMatrix", &CalcKINGMatrix, "");
  m.def("CalcUGRMMatrix", &CalcUGRMMatrix, "");
  m.def("FindUnrelatedGroup", &FindUnrelatedGroup, "");
  m.def("CalcUniLMM", &CalcUniLMM, "");
  m.def("CalcMultiLMM_REML", &CalcMultiLMM_REML, "");
  m.def("CalcMultiLMM_RML", &CalcMultiLMM_RML, "");
  m.def("CalcLinearRegressionGWAS", &CalcLinearRegressionGWAS, "");
  m.def("CalcLogisticGWAS", &CalcLogisticGWAS, "");
  m.def("CalcCCAGWAS", &CalcCCAGWAS, "");
  m.def("CalcUniLMMGWAS", &CalcUniLMMGWAS, "");
  m.def("UpdateAf", &UpdateAf, "");
  m.def("GenerateIndividuals", &GenerateIndividuals, "");
  m.def("GenerateAdmixedIndividuals", &GenerateAdmixedIndividuals, "");
  m.def("GeneratePairwiseSiblings", &GeneratePairwiseSiblings, "");
}