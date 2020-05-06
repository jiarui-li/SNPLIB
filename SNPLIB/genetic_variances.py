import math
import re
import numpy as np
import numpy.linalg as npl
from scipy import stats
import _SNPLIB as lib


def CalcUniLMM(traits, covariates, relationship_matrix, num_threads):
    w, v = npl.eig(relationship_matrix)
    std = np.std(traits, axis=0, dtype='double', ddof=1)
    Y = v.T@stats.zscore(traits, axis=0, ddof=1)
    X = v.T@covariates
    Vars, Res = lib.CalcUniLMM(Y, X, w, num_threads)
    Vars = Vars*(std**2)
    Res = Res*std
    return Vars, Res


def CalcMultiLMM_REML(traits, covariates, relationship_matrix, num_dims, num_threads):
    w, v = npl.eig(relationship_matrix)
    std = np.std(traits, axis=0, dtype='double', ddof=1)
    Y = v.T@stats.zscore(traits, axis=0, ddof=1)
    X = v.T@covariates
    Vars, Res = lib.CalcUniLMM(Y, X, w, num_threads)
    sigma_e, sigma_g = lib.CalcMultiLMM_REML(
        Y, X, w, Res, Vars, num_dims, num_threads)
    num_traits = Y.shape[1]//num_dims
    for i in range(num_traits):
        D = np.diag(std[i*num_dims:(i+1)*num_dims])
        sigma_e[..., i] = D@sigma_e[..., i]@D
        sigma_g[..., i] = D@sigma_g[..., i]@D
    return sigma_e, sigma_g


def CalcMultiLMM_RML(traits, covariates, relationship_matrix, num_dims, num_threads):
    w, v = npl.eig(relationship_matrix)
    std = np.std(traits, axis=0, dtype='double', ddof=1)
    Y = v.T@stats.zscore(traits, axis=0, ddof=1)
    X = v.T@covariates
    Vars, Res = lib.CalcUniLMM(Y, X, w, num_threads)
    sigma_e, sigma_g = lib.CalcMultiLMM_RML(
        w, Res, Vars, num_dims, num_threads)
    num_traits = Y.shape[1]//num_dims
    for i in range(num_traits):
        D = np.diag(std[i*num_dims:(i+1)*num_dims])
        sigma_e[..., i] = D@sigma_e[..., i]@D
        sigma_g[..., i] = D@sigma_g[..., i]@D
    return sigma_e, sigma_g
