function loadings = CalcPCALoadingsApprox(obj, nComponents, nParts)
if nargin<3
    nParts = 10;
end
L = 2*nComponents;
I = 10;
G = normrnd(0,1,[obj.nSamples,L]);
G = [G,zeros(obj.nSamples, I*L)];
af = obj.CalcAlleleFrequency();
grm = obj.CalcGRMMatrix();
for i=1:I
    G(:,i*L+1:(i+1)*L) = grm*G(:,(i-1)*L+1:i*L);
end
H = zeros(obj.nSNPs, L*(I+1));
nSNPsParts = ceil(obj.nSNPs/nParts);
for i=1:nParts-1
    ind = (i-1)*nSNPsParts+1:i*nSNPsParts;
    A = UnpackPCA_(obj.GENO(:,ind),af(ind),obj.nSamples)';
    H(ind,:) = A*G;
end
ind = (nParts-1)*nSNPsParts+1:obj.nSNPs;
A = UnpackPCA_(obj.GENO(:,ind),af(ind),obj.nSamples)';
H(ind,:) = A*G;
[Q,~] = qr(H,0);
T = zeros(obj.nSamples,L*(I+1));
for i=1:nParts-1
    ind = (i-1)*nSNPsParts+1:i*nSNPsParts;
    A = UnpackPCA_(obj.GENO(:,ind),af(ind),obj.nSamples);
    T = T + A*Q(ind,:);
end
ind = (nParts-1)*nSNPsParts+1:obj.nSNPs;
A = UnpackPCA_(obj.GENO(:,ind),af(ind),obj.nSamples);
T = T + A*Q(ind,:);
[~,S,W] = svd(T,'econ');
U = Q*W;
S = S(1:nComponents,1:nComponents);
U = U(:,1:nComponents);
loadings = S\U';
end

