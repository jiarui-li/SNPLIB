function [betas,p_values] = CalcUniMLMGWAS(obj,trait,covariates,relationship_matrix,block_size)
%CALCUNIMLMGWAS Summary of this function goes here
%   Detailed explanation goes here
if nargin<5
    block_size = 12000 * obj.nThreads;
end
if ~issymmetric(relationship_matrix)
    error("Please provide a symmetric relationship matrix");
    return
end
af = obj.CalcAlleleFrequency();
[V,D] = eig(relationship_matrix,'vector');
[~,ind] = sort(D,'descend');
D = D(ind);
V = V(:,ind);
X = V'*covariates;
Y = V'*zscore(trait);
num_blocks = ceil(obj.nSNPs / block_size);
p_values = zeros(obj.nSNPs,1);
betas = zeros(obj.nSNPs,1);
for i=1:num_blocks-1
    ind = (i-1)*block_size+1:i*block_size;
    geno_d = UnpackGENO_(obj.GENO(:,ind),obj.nSamples,af(ind));
    geno_d = V'*geno_d;
    [b,f,d] = CalcUniMLMGWAS_(Y,geno_d,X,D,obj.nThreads);
    p = fcdf(f,1,d,'upper');
    betas(ind) = b;
    p_values(ind) = p;
end
ind = (num_blocks-1)*block_size+1:obj.nSNPs;
geno_d = UnpackGENO_(obj.GENO(:,ind),obj.nSamples,af(ind));
[b,f,d] = CalcUniMLMGWAS_(Y,geno_d,X,D,obj.nThreads);
p = fcdf(f,1,d,'upper');
betas(ind) = b;
p_values(ind) = p;
end

