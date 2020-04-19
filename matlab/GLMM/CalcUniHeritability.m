function [h2_f, h2_r] = CalcUniHeritability(traits,covariates,genomic_background,relationship_matrix,num_threads)
%CALCUNIHERITABILITY Summary of this function goes here
%   Detailed explanation goes here
if nargin<5
    num_threads = feature('numcores');
end
if ~issymmetric(relationship_matrix)
    error("Please provide a symmetric relationship matrix");
    return
end
num_samples = size(traits,1);
[V,D] = eig(relationship_matrix,'vector');
[~,ind] = sort(D,'descend');
D = D(ind);
V = V(:,ind);
cov = [ones(num_samples,1),covariates,genomic_background];
X = V'*cov;
Y = V'*zscore(traits);
[vars, ~] = CalcUniMLM_(Y,X,D,num_threads);
h2_r = vars(2,:).^2;
Y = zscore(traits);
X = cov;
Hat = eye(size(Y,1))-X/(X'*X)*X';
res = Hat*Y;
res = res.^2;
v = sum(res,1);
X = [ones(num_samples,1),covariates];
Hat = eye(size(Y,1))-X/(X'*X)*X';
res = Hat*Y;
res = res.^2;
v = sum(res,1);
h2_rss_f = v/(num_samples-size(X,2))-vars(1,:).^2;
end