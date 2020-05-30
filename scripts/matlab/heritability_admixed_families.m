clear all;
%% Simulation settings
Fst = 0.2;
num_generations = 200;
effective_sample_size = floor(num_generations/2/(1-exp(-Fst)));
num_snps = 100000;
num_causal_snps = 10000;
num_samples = 4000;
num_traits = 2000;
%% Simulating Population Structure
aaf = unifrnd(0.1,0.9,[1,num_snps]);
pop_af = UpdateAf_(aaf,2,num_generations,effective_sample_size)';
%pop_af = [pop_af,repmat(aaf(num_snps/2+1:end),[2 1])];
pop_af(:,1) = [0.1;0.9];
pop = zeros(num_samples, 2);
pop(1:num_samples/2,1) = betarnd(0.1,0.1,[num_samples/2 1]);
pop(1:num_samples/2,2) = 1- pop(1:num_samples/2,1);
for i=1:num_samples/2
    pop(num_samples/2+i,1) = (pop(2*(i-1)+1,1)+pop(2*i,1))/2;
    pop(num_samples/2+i,2) = 1-pop(num_samples/2+i,1);
end
af = pop(1:num_samples/2,:)*pop_af;
%% Simulating Genotypes
obj = SNPLIB();
obj.nThreads = 4;
geno_p = GenerateIndividuals_(af);
geno_s = GeneratePairwiseSiblings_(geno_p,num_samples/2);
GENO = [geno_p;geno_s];
obj.CreateFromSIM(GENO,num_samples);
true_grm = obj.CalcAdmixedGRMMatrix(pop_af,pop);
grm = obj.CalcGRMMatrix();
scores = obj.CalcSUGIBSScores(1);
%% Simulating traits
true_genetic_corr = zeros(num_pairs,1);
true_env_corr = zeros(num_pairs,1);
sim_traits = zeros(num_samples, num_traits);
all_snps = 2:num_snps;
snp_ind = randsample(all_snps,num_causal_snps,false);
obj.extract([1,snp_ind]);
geno_d = obj.UnpackGeno();
parfor i=1:num_traits
    random_effect = geno_d(:,2:end)*normrnd(0,1,[num_causal_snps 1]);
    env_effect = zscore(normrnd(0,1,[num_samples,1]));
    gene_effect = zscore(random_effect)*sqrt(0.5)+zscore(geno_d(:,1))*sqrt(0.5);
    env_effect = env_effect-dot(env_effect,gene_effect)/dot(gene_effect,gene_effect).*gene_effect;
    sim_traits(:,i) = zscore(gene_effect)*sqrt(0.72)+zscore(env_effect)*sqrt(0.28); 
end
[~,h2_r_1] = CalcUniHeritability(sim_traits,[],[],grm,8);
[h2_f_2,h2_r_2] = CalcUniHeritability(sim_traits,[],scores,grm,8);
[h2_f,h2_r] = CalcUniHeritability(sim_traits,[],scores,true_grm,8);