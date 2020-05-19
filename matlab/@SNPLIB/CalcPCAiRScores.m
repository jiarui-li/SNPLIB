function scores = CalcPCAiRScores(obj,nComp,nParts)
%CALCPCAIRSCORES Calculate PCA scores with relatives
%   Detailed explanation goes here
if nargin < 3
    nParts = 10;
end
ind_u = obj.FindUnrelatedGroup()';
geno_u = Keep_(obj.GENO,obj.nSamples,int32(ind_u-1));
nSamples_u = length(ind_u);
af = CalcAlleleFrequencies_(geno_u,nSamples_u);
A = UnpackGRMGeno_(geno_u,af,nSamples_u)';
[U,S,~] = svds(A,nComp);
loadings = S\U';
scores = zeros(nComp,obj.nSamples);
nSNPsParts = ceil(obj.nSNPs/nParts);
for i=1:nParts-1
    ind = (i-1)*nSNPsParts+1:i*nSNPsParts;
    A = UnpackGRMGeno_(obj.GENO(:,ind),af(ind),obj.nSamples)';
    scores = scores + loadings(:,ind)*A;
end
ind = (nParts-1)*nSNPsParts+1:obj.nSNPs;
A = UnpackGRMGeno_(obj.GENO(:,ind),af(ind),obj.nSamples)';
scores = scores + loadings(:,ind)*A;
scores = scores';
end

