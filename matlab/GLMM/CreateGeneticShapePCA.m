function out = CreateGeneticShapePCA(landmarks,covariates,relationship_matrix,num_components,num_threads)
if nargin<5, num_threads = feature('numcores'); end
if nargin<4, num_components = 0; end
if ~issymmetric(relationship_matrix)
    error("Please provide a symmetric relationship matrix");
    return
end

[V,D] = eig(relationship_matrix,'vector');
[~,ind] = sort(D,'descend');
D = D(ind);
V = V(:,ind);
cov = [ones(size(traits,1),1),covariates];
out.rm_eigenvectors = V;
out.rm_eigenvalues = D;
s = std(landmarks);
X = V'*cov;
Y = V'*zscore(landmarks);
disp('Calculating Shape PCA');
[~,res] = CalcUniMLM_(Y,X,D,num_threads);
res = V*res;
res = bxsfun(@times,res,s);
num_dims = size(landmarks,2);
[~,S,loadings] = svd(res,'econ');
S = diag(S).^2;
if num_components == 0
    s_shuffle = zeros(length(S), 100);
    res_shuffle = res;
    for i=1:100
        for j=1:num_dims
            res_shuffle(:,j) = res(randperm(num_samples),j);
        end
        s_shuffle(:,i) = svd(res_shuffle);
    end
    s_shuffle = s_shuffle.^2;
    s_low = quantile(s_shuffle',0.05);
    num_components = find(S<s_low',1,'first');
end
out.loadings = loadings(:,1:num_components);
C = cumsum(S)/sum(S);
out.var_explained = C(num_components);
out.num_components;
disp([num2str(num_components),' ',num2str(var_explained)]);

end