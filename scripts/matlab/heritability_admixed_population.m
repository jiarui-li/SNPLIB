%% Simulation settings
Fst = 0.1;
num_generations = 200;
effective_sample_size = floor(num_generations/2/(1-exp(-Fst)));
num_snps = 100000;
num_causal_snps = 10000;
num_samples = 4000;
num_traits = 5000;
num_pairs = num_traits/2;
%% Simulating Population Structure
aaf = unifrnd(0.05,0.95,[1,num_snps]);
pop_af = UpdateAf_(aaf,2,num_generations,effective_sample_size);
pop = zeros(num_samples, 2);
pop(:,1) = sort(betarnd(0.1,0.1,[num_samples 1]));
pop(:,2) = 1- pop(:,1);
af = pop*pop_af;
%% Simulating Genotypes
obj = SNPLIB();
obj.nThreads = 4;
obj.GenerateIndividuals(af);
geno_d = obj.UnpackGeno();
%% Simulating traits
true_genetic_corr = zeros(num_pairs,1);
true_env_corr = zeros(num_pairs,1);
sim_traits = zeros(num_samples, num_traits);
parfor i=1:num_pairs
    all_snps = 1:num_snps;
    additive_snp = randsample(all_snps,1);
    all_snps = setdiff(all_snps,additive_snp);
    snp_ind1 = randsample(all_snps,num_causal_snps,false);
    all_snps = setdiff(all_snps,snp_ind1);
    num_shared_snps = floor(rand*num_causal_snps);
    snp_ind2 = randsample(snp_ind1,num_shared_snps,false);
    snp_ind2 = [snp_ind2,randsample(all_snps,num_causal_snps-num_shared_snps,false)];
    rho = rand*2-1;
    env_effects = mvnrnd([0,0],[1,rho;rho,1],num_samples);
    
end