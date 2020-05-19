function [Sigma_e,Sigma_g] = CalcMultiLMM_REML(traits,covariates,relationship_matrix,num_dims,num_threads)
if nargin<5
    num_threads = feature('numcores');
end
[V,D] = eig(relationship_matrix,'vector');
[~,ind] = sort(D,'descend');
D = D(ind);
V = V(:,ind);
cov = [ones(size(traits,1),1),covariates,genomic_background];
X = V'*cov;
Y = V'*zscore(traits);
[vars,res] = CalcUniMLM_(Y,X,D,num_threads);
[Sigma_e,Sigma_g] = CalcMultiLMM_REML_(Y,X,D,res,vars,num_dims,num_threads);
end

