function [h2_f,h2_r] = CalcUniHeritability(traits,covariates,genomic_background,relationship_matrix,num_threads)
if nargin<5
    num_threads = feature('numcores');
end
num_samples = size(traits,1);
[V,D] = eig(relationship_matrix,'vector');
[~,ind] = sort(D,'descend');
D = D(ind);
V = V(:,ind);
cov = [ones(num_samples,1),covariates,genomic_background];
X = V'*cov;
Y = V'*zscore(traits);
[vars,~] = CalcUniLMM_(Y,X,D,num_threads);
h2_r = vars(2,:).^2;
Y = zscore(traits);
X = [ones(num_samples,1),covariates];
Hat = eye(size(Y,1))-X/(X'*X)*X';
res = Hat*Y;
res = res.^2;
v = sum(res,1);
h2_f = v/(num_samples-size(X,2))-vars(1,:).^2;
end

