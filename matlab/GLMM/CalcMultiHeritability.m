function [h2_f,h2_r] = CalcMultiHeritability(traits,covariates,genomic_background,relationship_matrix,num_dims,num_threads)
%CALCMULTIHERITABILITY Summary of this function goes here
%   Detailed explanation goes here
if nargin<6
    num_threads = feature('numcores');
end
if ~issymmetric(relationship_matrix)
    error("Please provide a symmetric relationship matrix");
    return
end
Y = traits;
X = [ones(size(traits,1),1),covariates];
num_cov = size(X,2);
Hat = eye(size(Y,1))-X/(X'*X)*X';
res = Hat*Y;
res = res.^2;
v = sum(res)./(size(Y,1)-num_cov);
[V,D] = eig(relationship_matrix,'vector');
[~,ind] = sort(D,'descend');
D = D(ind);
V = V(:,ind);
cov = [ones(size(traits,1),1),covariates,genomic_background];
X = V'*cov;
Y = V'*zscore(traits);
[vars,res] = CalcUniMLM_(Y,X,D,num_threads);
num_traits = size(traits,2);
num_traits = num_traits / num_dims;
[Sigma_e,Sigma_g] = CalcMMLMSigmas_(Y,X,D,vars,res,num_dims,num_threads);
h2_r = zeros(num_traits,1);
h2_f = zeros(num_traits,1);
s = std(traits);
for i=1:num_traits
    ind = (i-1)*num_dims+1:i*num_dims;
    S = s(ind)'*s(ind);
    v_e = S.*Sigma_e(:,:,i);
    v_g = S.*Sigma_g(:,:,i);
    h2_r(i) = sum(diag(v_g))/sum(s(ind).^2);
    h2_f(i) = (sum(v(ind))-sum(diag(v_e)))/sum(s(ind).^2);
end
end

