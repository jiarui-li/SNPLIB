function scores = ProjectUPCA(obj,loadings,nParts)
if nargin < 3
    nParts = 10;
end
nComp = size(loadings,1);
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
scores = scores';
end