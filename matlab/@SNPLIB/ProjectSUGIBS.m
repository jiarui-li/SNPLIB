function scores = ProjectSUGIBS(obj,ref_obj,loadings,nParts)
if nargin < 4
    nParts = 10;
end
if ~isa(ref_obj,'SNPLIB')
    error('Please provide a SNPLIB class object as reference!');
end
nComp = size(loadings,1);
scores = zeros(nComp,obj.nSamples);
nSNPsParts = ceil(obj.nSNPs/nParts);
for i=1:nParts-1
    ind = (i-1)*nSNPsParts+1:i*nSNPsParts;
    A = UnpackUG_(obj.GENO(:,ind),obj.nSamples)';
    scores = scores + loadings(:,ind)*A;
end
ind = (nParts-1)*nSNPsParts+1:obj.nSNPs;
A = UnpackUG_(obj.GENO(:,ind),obj.nSamples)';
scores = scores + loadings(:,ind)*A;
connect = CalcIBSConnect_(ref_obj.GENO,ref_obj.nSamples,obj.GENO,obj.nSamples,obj.nThreads);
scores = scores*diag(connect.^(-1));
scores = scores';
end