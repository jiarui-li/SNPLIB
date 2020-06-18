"""SNPLIB class Module

This module contains the definition of SNPLIB class and some supportive methods.

Example
-------

Notes
-----
"""
import math
import re
import numpy as np
import numpy.linalg as npl
import pandas as pd
from copy import deepcopy
from scipy.stats import t, chi2, f
from multiprocessing import cpu_count
import _SNPLIB as lib


chr_dict = {
    'X': 23,
    'Y': 24,
    'XY': 25,
    'MT': 26,
}


def convert_chr(chr_c):
    if chr_c.isdigit():
        return int(chr_c)
    else:
        return chr_dict[chr_c.upper()]


def CalcIBSConnection(src_geno, dest_geno, num_threads):
    return lib.CalcIBSConnection(src_geno.GENO, dest_geno.GENO, src_geno.nSamples, dest_geno.nSamples, num_threads)


def UpdateAf(aaf, num_pops, num_generations, effective_sample_size):
    return lib.UpdateAf(aaf, num_pops, num_generations, effective_sample_size).T


class SNPLIB:
    """The SNPLIB class stores the genotype data in plink binary format and provides methods analyzing genotype data

    Attributes
    ----------
    nThreads : int
        Number of threads used for computation
    SNPs : :obj:`pandas.DataFrame`, optional
        `pandas.DataFrame` of SNPs' information imported from plink .bim file
    Samples : :obj:`pandas.DataFrame`, optional
        `pandas.DataFrame` of Samples' information imported from plink .fam file
    """

    def __init__(self, nThreads=cpu_count()//2):
        """Constructor function

        Parameters
        ----------
        nThreads : int, optional
            Number of threads used for computation

        """
        self.nThreads = nThreads
        self.SNPs = []
        self.Samples = []

    def importPLINKDATA(self, bfile):
        """Import plink binary fileset

        Parameters
        ----------
        bfile : str
            The name of plink binary fileset

        """
        filename = bfile + '.bim'
        self.SNPs = pd.read_table(
            bfile+'.bim', sep=None, names=['CHR', 'RSID', 'Cm', 'POS', 'ALT', 'REF'], engine='python')
        self.Samples = pd.read_table(bfile+'.fam', sep=None,
                                     names=['FID', 'IID', 'PID', 'MID', 'Sex', 'Pheno'], engine='python')
        self.nSNPs = self.SNPs.shape[0]
        self.nSamples = self.Samples.shape[0]
        filename = bfile + '.bed'
        num_bytes = math.ceil(self.nSamples / 4.0)
        GENO = np.fromfile(filename, dtype=np.uint8, count=-1)
        GENO = GENO[3:]
        self.GENO = np.reshape(GENO, (num_bytes, - 1), order='F')

    # Simulations
    def GenerateIndividuals(self, af):
        """Simulate the genotypes according to the individual allele frequencies

        Parameters
        ----------
        af : ndarray
            A `ndarray` matrix contains the individual allele frequencies, with shape of  ``(num_samples, num_snps)``
        """
        self.nSamples = af.shape[0]
        self.nSNPs = af.shape[1]
        self.GENO = lib.GenerateIndividuals(af)

    def GenerateAdmixedIndividuals(self, af, num_samples):
        self.nSamples = num_samples
        self.nSNPs = af.shape[1]
        self.GENO = lib.GenerateAdmixedIndividuals(af, num_samples)

    def GeneratePairwiseSiblings(self, parent_obj):
        self.nSamples = parent_obj.nSamples
        self.nSNPs = parent_obj.nSNPs
        self.GENO = lib.GeneratePairwiseSiblings(
            parent_obj.GENO, self.nSamples//2)

    # Data manage
    def UnpackGeno(self):
        return lib.UnpackGeno(self.GENO, self.nSamples)

    def Keep(self, index):
        result = SNPLIB(self.nThreads)
        result.nSamples = len(index)
        result.nSNPs = deepcopy(self.nSNPs)
        result.GENO = lib.Keep(self.GENO, self.nSamples, index)
        return result

    def Extract(self, index):
        result = SNPLIB(self.nThreads)
        result.nSamples = deepcopy(self.nSamples)
        result.nSNPs = len(index)
        result.GENO = deepcopy(self.GENO[:, index])
        return result
    # Statistics

    def CalcAlleleFrequencies(self):
        return lib.CalcAlleleFrequencies(self.GENO, self.nSamples)

    def CalcMissing(self):
        return lib.CalcMissing(self.GENO, self.nSamples)

    def CalcAdjustedAF(self, covariates):
        return lib.CalcAdjustedAF(self.GENO, covariates, self.nThreads)

    def CalcAdjustedMAF(self, covariates):
        return lib.CalcAdjustedMAF(self.GENO, covariates, self.nThreads)

    # Relationships
    def CalcAdjustedGRM(self, covariates):
        matrix, gcta_diag = lib.CalcAdjustedGRM(
            self.GENO, covariates, self.nThreads)
        return matrix, gcta_diag

    def CalcAdmixedGRM(self, pop_af, pop):
        matrix, gcta_diag = lib.CalcAdmixedGRM(
            self.GENO, pop_af, pop, self.nThreads)
        return matrix, gcta_diag

    def CalcGRMMatrix(self):
        af = self.CalcAlleleFrequencies()
        matrix = lib.CalcGRMMatrix(self.GENO, af, self.nSamples, self.nThreads)
        gcta_diag = lib.CalcGCTADiagonal(
            self.GENO, af, self.nSamples, self.nThreads)
        return matrix, gcta_diag

    def CalcIBSMatrix(self):
        return lib.CalcIBSMatrix(self.GENO, self.nSamples, self.nThreads)

    def CalcKINGMatrix(self):
        return lib.CalcKINGMatrix(self.GENO, self.nSamples, self.nThreads)

    def CalcUGRMMatrix(self):
        return lib.CalcUGRMMatrix(self.GENO, self.nSamples, self.nThreads)

    def FindUnrelated(self, threshold=0.044):
        king = self.CalcKINGMatrix()
        return lib.FindUnrelatedGroup(king, threshold)

    # Ancestry Estimation

    def CalcPCAScores(self, nComponents):
        grm, _ = self.CalcGRMMatrix()
        w, V = npl.eig(grm)
        ind = np.argsort(-w)
        return V[:, ind[:nComponents]]

    def CalcPCALoadingsExact(self, nComponents):
        af = self.CalcAlleleFrequencies()
        A = lib.UnpackGRMGeno(self.GENO, af, self.nSamples)
        U, s, _ = npl.svd(A.T, full_matrices=False)
        S = np.diag(s[:nComponents])
        U = U[:, :nComponents]
        return npl.solve(S, U.T)

    def CalcPCALoadingsApprox(self, nComponents, nParts=10):
        af = self.CalcAlleleFrequencies()
        grm = lib.CalcGRMMatrix(self.GENO, af, self.nSamples, self.nThreads)
        L = 2*nComponents
        I = 10
        G = np.zeros((self.nSamples, L*(I+1)), dtype='double', order='F')
        G[:, :L] = np.random.randn(self.nSamples, L)
        for i in range(1, I):
            G[:, i*L: (i+1)*L] = grm@G[:, (i-1)*L:i*L]

        H = np.zeros((self.nSNPs, L*(I+1)), dtype='double', order='F')
        nSNPsPart = math.ceil(self.nSNPs/nParts)
        for i in range(nParts-1):
            A = lib.UnpackGRMGeno(
                self.GENO[:, i*nSNPsPart:(i+1)*nSNPsPart], af[i*nSNPsPart:(i+1)*nSNPsPart], self.nSamples)
            H[i*nSNPsPart:(i+1)*nSNPsPart, :] = A.T@G

        A = lib.UnpackGRMGeno(
            self.GENO[:, (nParts-1)*nSNPsPart:], af[(nParts-1)*nSNPsPart:], self.nSamples)
        H[(nParts-1)*nSNPsPart:, :] = A.T@G
        Q, _ = npl.qr(H)
        T = np.zeros((self.nSamples, L*(I+1)), dtype='double', order='F')
        for i in range(nParts-1):
            A = lib.UnpackGRMGeno(
                self.GENO[:, i*nSNPsPart:(i+1)*nSNPsPart], af[i*nSNPsPart:(i+1)*nSNPsPart], self.nSamples)
            T = T + A@Q[i*nSNPsPart:(i+1)*nSNPsPart, :]

        A = lib.UnpackGRMGeno(
            self.GENO[:, (nParts-1)*nSNPsPart:], af[(nParts-1)*nSNPsPart:], self.nSamples)
        T = T+A@Q[(nParts-1)*nSNPsPart:, :]
        _, S, W = npl.svd(T, full_matrices=False)
        U = Q@W.T
        S = S[:nComponents]
        S = np.diag(S)
        U = U[:, :nComponents]
        return npl.solve(S, U.T)

    def ProjectPCA(self, ref_obj, loadings, nParts=10):
        nSNPsPart = math.ceil(self.nSNPs/nParts)
        af = ref_obj.CalcAlleleFrequencies()
        nComponents = loadings.shape[0]
        scores = np.zeros((nComponents, self.nSamples))
        for i in range(nParts-1):
            A = lib.UnpackGRMGeno(
                self.GENO[:, i*nSNPsPart:(i+1)*nSNPsPart], af[i*nSNPsPart:(i+1)*nSNPsPart], self.nSamples)
            scores = scores + loadings[:, i*nSNPsPart:(i+1)*nSNPsPart]@A.T

        A = lib.UnpackGRMGeno(
            self.GENO[:, (nParts-1)*nSNPsPart:], af[nParts*nSNPsPart:], self.nSamples)
        scores = scores + loadings[:, (nParts-1)*nSNPsPart:]@A.T
        return scores.T

    def CalcPCAiRScores(self, nComponents, is_exact=False):
        ind = self.FindUnrelated()
        Unrelated = self.Keep(ind)
        if is_exact:
            loadings = Unrelated.CalcPCALoadingsExact(nComponents)
        else:
            loadings = Unrelated.CalcPCALoadingsApprox(nComponents)
        return self.ProjectPCA(Unrelated, loadings)

    def CalcSUGIBSScores(self, nComponents):
        ibs = self.CalcIBSMatrix()
        d = np.sum(ibs, axis=0)
        ugrm = self.CalcUGRMMatrix()
        D = np.diag(d**-0.5)
        I = D@ugrm@D
        w, V = npl.eig(I)
        ind = np.argsort(-w)
        return D@V[:, ind[1:nComponents+1]]

    def CalcSUGIBSLoadingsExact(self, nComponents):
        ibs = self.CalcIBSMatrix()
        d = np.sum(ibs, axis=0)
        D = np.diag(d**-0.5)
        A = lib.UnpackUGeno(self.GENO, self.nSamples)
        A = D@A
        U, s, _ = npl.svd(A.T, full_matrices=False)
        S = np.diag(s[1:nComponents+1])
        U = U[:, 1:nComponents+1]
        return npl.solve(S, U.T)

    def CalcSUGIBSLoadingsApprox(self, nComponents, nParts=10):
        L = 2*(nComponents+1)
        I = 10
        G = np.zeros((self.nSamples, L*(I+1)), dtype='double', order='F')
        G[:, :L] = np.random.randn(self.nSamples, L)
        ibs = self.CalcIBSMatrix()
        ugrm = self.CalcUGRMMatrix()
        d = np.sum(ibs, axis=0)
        D = np.diag(d**-0.5)
        ugrm = D@ugrm@D
        for i in range(1, I):
            G[:, i*L: (i+1)*L] = ugrm@G[:, (i-1)*L:i*L]

        H = np.zeros((self.nSNPs, L*(I+1)), dtype='double', order='F')
        nSNPsPart = math.ceil(self.nSNPs/nParts)
        for i in range(nParts-1):
            A = lib.UnpackUGeno(
                self.GENO[:, i*nSNPsPart:(i+1)*nSNPsPart], self.nSamples)
            H[i*nSNPsPart:(i+1)*nSNPsPart, :] = A.T@D@G

        A = lib.UnpackUGeno(self.GENO[:, (nParts-1)*nSNPsPart:], self.nSamples)
        H[(nParts-1)*nSNPsPart:, :] = A.T@D@G
        Q, _ = npl.qr(H)
        T = np.zeros((self.nSamples, L*(I+1)), dtype='double', order='F')
        for i in range(nParts-1):
            A = lib.UnpackUGeno(
                self.GENO[:, i*nSNPsPart:(i+1)*nSNPsPart], self.nSamples)
            T = T + D@A@Q[i*nSNPsPart:(i+1)*nSNPsPart, :]

        A = lib.UnpackUGeno(self.GENO[:, (nParts-1)*nSNPsPart:], self.nSamples)
        T = T+D@A@Q[(nParts-1)*nSNPsPart:, :]
        _, S, W = npl.svd(T, full_matrices=False)
        U = Q@W.T
        S = S[1:nComponents+1]
        S = np.diag(S)
        U = U[:, 1:nComponents+1]
        return npl.solve(S, U.T)

    def ProjectSUGIBS(self, ref_obj, loadings, nParts=10):
        nSNPsPart = math.ceil(self.nSNPs/nParts)
        nComponents = loadings.shape[0]
        scores = np.zeros((nComponents, self.nSamples))
        for i in range(nParts-1):
            A = lib.UnpackUGeno(
                self.GENO[:, i*nSNPsPart:(i+1)*nSNPsPart], self.nSamples)
            scores = scores + loadings[:, i*nSNPsPart:(i+1)*nSNPsPart]@A.T

        A = lib.UnpackUGeno(
            self.GENO[:, (nParts-1)*nSNPsPart:], self.nSamples)
        scores = scores + loadings[:, (nParts-1)*nSNPsPart:]@A.T
        connect = CalcIBSConnection(ref_obj, self, self.nThreads)
        D = np.diag(connect**-1)
        scores = scores@D
        return scores.T

    def CalcSUGIBSiRScores(self, nComponents, is_exact=False):
        ind = self.FindUnrelated()
        Unrelated = self.Keep(ind)
        if is_exact:
            loadings = Unrelated.CalcSUGIBSLoadingsExact(nComponents)
        else:
            loadings = Unrelated.CalcSUGIBSLoadingsApprox(nComponents)
        return self.ProjectSUGIBS(Unrelated, loadings)

        # GWAS
    def CalcLinearRegressionGWAS(self, trait, covariates):
        betas, stats = lib.CalcLinearRegressionGWAS(
            self.GENO, covariates, trait, self.nThreads)
        df = len(trait)-covariates.shape[1]
        pvalues = t.cdf(stats, df=df)
        return betas, pvalues

    def CalcCCAGWAS(self, trait):
        Y = trait - trait.mean(axis=0)
        betas, rho2 = lib.CalcCCAGWAS(self.GENO, Y, self.nThreads)
        nDims = trait.shape[1]
        t = (nDims**2-4)/(nDims**2+1-5)
        Lambda = 1.0 - rho2
        t = np.sqrt(t)
        w = self.nSamples - (nDims+4)/2
        df1 = nDims
        df2 = w*t-nDims/2+1
        Lambda = np.power(Lambda, 1.0/t)
        F = (1-Lambda)/Lambda*df2/df1
        pvalues = f.cdf(F, dfn=df2, dfd=df1)
        return betas, rho2, pvalues
