function scores = CalcSUGIBSiRScores(obj,nComp, nParts)
%CALCSUGIBSIRSCORES Calculate SUGIBS with relatives
%   Detailed explanation goes here
if nargin < 3
    nParts = 10;
end
ind_u = obj.FindUnrelated();
geno_u = Keep_(obj.GENO,obj.nSamples,int32(ind_u-1));
nSamples_u = length(ind_u);
ibs = CalcIBSMatrix_(geno_u,nSamples_u,obj.nThreads);
d = sum(ibs);
d = diag(d.^(-0.5));
A = UnpackUGeno_(geno_u,nSamples_u);
A = A'*d;
[U,S,~] = svds(A,nComp+1);
S = S(2:nComp+1,2:nComp+1);
U = U(:,2:nComp+1);
loadings = S\U';
scores = zeros(nComp,obj.nSamples);
nSNPsParts = ceil(obj.nSNPs/nParts);
for i=1:nParts-1
    ind = (i-1)*nSNPsParts+1:i*nSNPsParts;
    A = UnpackUGeno_(obj.GENO(:,ind),obj.nSamples)';
    scores = scores + loadings(:,ind)*A;
end
ind = (nParts-1)*nSNPsParts+1:obj.nSNPs;
A = UnpackUGeno_(obj.GENO(:,ind),obj.nSamples)';
scores = scores + loadings(:,ind)*A;
connect = CalcIBSConnection_(geno_u,nSamples_u,obj.GENO,obj.nSamples,obj.nThreads);
scores = scores*diag(connect.^(-1));
scores = scores';
end

