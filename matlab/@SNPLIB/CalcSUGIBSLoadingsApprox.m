function loadings = CalcSUGIBSLoadingsApprox(obj, nComponents, nParts)
if nargin<3
    nParts = 10;
end
L = 2*(nComponents+1);
I = 10;
G = normrnd(0,1,[obj.nSamples,L]);
G = [G,zeros(obj.nSamples, I*L)];
ibs = obj.CalcIBSMatrix();
d = sum(ibs);
d = diag(d.^(-0.5));
ugrm = CalcUGRMMatrix_(obj.GENO,obj.nSamples,obj.nThreads);
IBS = d*ugrm*d;
for i=1:I
    G(:,i*L+1:(i+1)*L) = IBS*G(:,(i-1)*L+1:i*L);
end
H = zeros(obj.nSNPs, L*(I+1));
nSNPsParts = ceil(obj.nSNPs/nParts);
for i=1:nParts-1
    ind = (i-1)*nSNPsParts+1:i*nSNPsParts;
    A = UnpackUGeno_(obj.GENO(:,ind),obj.nSamples)';
    H(ind,:) = A*d*G;
end
ind = (nParts-1)*nSNPsParts+1:obj.nSNPs;
A = UnpackUGeno_(obj.GENO(:,ind),obj.nSamples)';
H(ind,:) = A*d*G;
[Q,~] = qr(H,0);
T = zeros(obj.nSamples,L*(I+1));
for i=1:nParts-1
    ind = (i-1)*nSNPsParts+1:i*nSNPsParts;
    A = UnpackUGeno_(obj.GENO(:,ind),obj.nSamples);
    T = T + d*A*Q(ind,:);
end
ind = (nParts-1)*nSNPsParts+1:obj.nSNPs;
A = UnpackUGeno_(obj.GENO(:,ind),obj.nSamples);
T = T + d*A*Q(ind,:);
[~,S,W] = svd(T,'econ');
U = Q*W;
S = S(2:nComponents+1,2:nComponents+1);
U = U(:,2:nComponents+1);
loadings = S\U';
end