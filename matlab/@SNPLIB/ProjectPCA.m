function scores = ProjectPCA(obj,ref_obj,loadings,nParts)
if nargin < 3
    nParts = 10;
end
if ~isa(ref_obj,'SNPLIB')
    error('Please provide a SNPLIB class object as reference!');
end
nComp = size(loadings,1);
scores = zeros(nComp,obj.nSamples);
af = ref_obj.CalcAlleleFrequency();
nSNPsParts = ceil(obj.nSNPs/nParts);
for i=1:nParts-1
    ind = (i-1)*nSNPsParts+1:i*nSNPsParts;
    A = UnpackPCA_(obj.GENO(:,ind),af(ind),obj.nSamples)';
    scores = scores + loadings(:,ind)*A;
end
ind = (nParts-1)*nSNPsParts+1:obj.nSNPs;
A = UnpackPCA_(obj.GENO(:,ind),af(ind),obj.nSamples)';
scores = scores + loadings(:,ind)*A;
scores = scores';
end