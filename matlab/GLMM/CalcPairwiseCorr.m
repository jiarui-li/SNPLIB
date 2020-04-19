function [genetic_corr,env_corr] = CalcPairwiseCorr(traits,covariates,relationship_matrix,num_dims,num_threads)

if nargin<6
    num_threads = feature('numcores');
end
if ~issymmetric(relationship_matrix)
    error("Please provide a symmetric relationship matrix");
    return
end
block_size = num_threads * 150;
num_traits = size(traits,2) / num_dims;
num_pairs = num_traits*(num_traits-1)/2;
num_blocks = ceil(num_pairs / block_size);
num_samples = size(traits,1);
[V,D] = eig(relationship_matrix,'vector');
[~,ind] = sort(D,'descend');
D = D(ind);
V = V(:,ind);
cov = [ones(size(traits,1),1),covariates];
num_covariates = size(cov,2);
X = V'*cov;
s = std(traits);
Y = V'*bsxfun(@rdivide,traits,s);
[vars,res] = CalcUniMLM_(Y,X,D,num_threads);
mkdir('GC_WORKING');
cd('GC_WORKING');
save('METAINFO','num_dims','num_traits','num_blocks');
ii = 1;
jj = 2;
tic
for l=1:num_blocks-1
    init_vars = zeros(2,2*num_dims*block_size);
    Res = zeros(num_samples, 2*num_dims*block_size);
    for i=1:block_size
        ind = [(ii-1)*num_dims+1:ii*num_dims,(jj-1)*num_dims+1:jj*num_dims];
        Res(:,(i-1)*2*num_dims+1:i*2*num_dims) = res(:,ind);
        init_vars(:,(i-1)*2*num_dims+1:i*2*num_dims) = vars(:,ind);
        jj = jj + 1;
        if jj>num_traits
            ii = ii + 1;
            jj = ii + 1;
        end
    end
    [Sigma_e,Sigma_g] = CalcRMLMSigmas_(Res,D,init_vars,num_covariates,2*num_dims,num_threads);
    in.Sigma_e = Sigma_e;
    in.Sigma_g = Sigma_g;
    save(['block_input_',num2str(l)],'in');
    tElapsed = toc / 60;
    tRemained = (num_blocks-l)*tElapsed/l;
    disp([num2str(l),'/', num2str(num_blocks), ' blocks finished. ',num2str(tElapsed),' minutes elapsed and ', num2str(tRemained),' minutes remain']);
end
num_pairs_left = num_pairs - (num_blocks-1)*block_size;
init_vars = zeros(2, 2*num_dims*num_pairs_left);
Res = zeros(num_samples, 2*num_dims*num_pairs_left);
for i=1:num_pairs_left
    ind = [(ii-1)*num_dims+1:ii*num_dims,(jj-1)*num_dims+1:jj*num_dims];
    Res(:,(i-1)*2*num_dims+1:i*2*num_dims) = res(:,ind);
    init_vars(:,(i-1)*2*num_dims+1:i*2*num_dims) = vars(:,ind);
    jj = jj + 1;
    if jj>num_traits
        ii = ii + 1;
        jj = ii + 1;
    end
end
[Sigma_e,Sigma_g] = CalcRMLMSigmas_(Res,D,init_vars,num_covariates,2*num_dims,num_threads);
in.Sigma_e = Sigma_e;
in.Sigma_g = Sigma_g;
save(['block_input_',num2str(num_blocks)],'in');
ii = 1;
jj = 2;
genetic_corr = zeros(num_traits,num_traits);
env_corr = zeros(num_traits,num_traits);
for l=1:num_blocks-1
    load(['block_input_',num2str(l)]);
    for i=1:block_size
        ind = [(ii-1)*num_dims+1:ii*num_dims,(jj-1)*num_dims+1:jj*num_dims];
        ss = s(ind);
        v_e = in.Sigma_e(:,:,i);
        v_g = in.Sigma_g(:,:,i);
        v_e = v_e + triu(v_e', 1);
        v_g = v_g + triu(v_g', 1);
        v_e = (ss'*ss).*v_e;
        v_g = (ss'*ss).*v_g;
        VE1 = sum(sum(v_e(1:num_dims,1:num_dims).^2));
        VE2 = sum(sum(v_e(num_dims+1:end,num_dims+1:end).^2));
        CE = sum(sum(v_e(1:num_dims,num_dims+1:end).^2));
        VG1 = sum(sum(v_g(1:num_dims,1:num_dims).^2));
        VG2 = sum(sum(v_g(num_dims+1:end,num_dims+1:end).^2));
        CG = sum(sum(v_g(1:num_dims,num_dims+1:end).^2));
        genetic_corr(ii,jj) = CG/sqrt(VG1*VG2);
        env_corr(ii,jj) = CE/sqrt(VE1*VE2);
        jj = jj + 1;
        if jj>num_traits
            ii = ii + 1;
            jj = ii + 1;
        end
    end
end
load(['block_input_',num2str(num_blocks)]);
for i=1:num_pairs_left
    ind = [(ii-1)*num_dims+1:ii*num_dims,(jj-1)*num_dims+1:jj*num_dims];
    ss = s(ind);
    v_e = in.Sigma_e(:,:,i);
    v_g = in.Sigma_g(:,:,i);
    v_e = v_e + triu(v_e', 1);
    v_g = v_g + triu(v_g', 1);
    v_e = (ss'*ss).*v_e;
    v_g = (ss'*ss).*v_g;
    VE1 = sum(sum(v_e(1:num_dims,1:num_dims).^2));
    VE2 = sum(sum(v_e(num_dims+1:end,num_dims+1:end).^2));
    CE = sum(sum(v_e(1:num_dims,num_dims+1:end).^2));
    VG1 = sum(sum(v_g(1:num_dims,1:num_dims).^2));
    VG2 = sum(sum(v_g(num_dims+1:end,num_dims+1:end).^2));
    CG = sum(sum(v_g(1:num_dims,num_dims+1:end).^2));
    genetic_corr(ii,jj) = CG/sqrt(VG1*VG2);
    env_corr(ii,jj) = CE/sqrt(VE1*VE2);
    jj = jj + 1;
    if jj>num_traits
        ii = ii + 1;
        jj = ii + 1;
    end
end
genetic_corr = genetic_corr + genetic_corr' + eye(num_traits);
env_corr = env_corr + env_corr' + eye(num_traits);
end