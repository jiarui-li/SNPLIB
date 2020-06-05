clear all;
%% Simulation settings
Fst = 0.2;
num_generations = 20;
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
pop(1:num_samples/2,1) = sort(betarnd(0.1,0.1,[num_samples/2 1]));
for i=1:num_samples/4
    pop(num_samples/2+2*(i-1)+1,1) = (pop(2*(i-1)+1,1)+pop(2*i,1))/2;
    pop(num_samples/2+2*i,1) = (pop(2*(i-1)+1,1)+pop(2*i,1))/2;
end
pop(:,2) = 1-pop(:,1);
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
scores = obj.CalcSUGIBSiRScores(1);
%% Simulating traits
sim_traits = zeros(num_samples, num_traits);
all_snps = 2:num_snps;
snp_ind = randsample(all_snps,num_causal_snps,false);
obj.extract([1,snp_ind]);
true_causal_grm = obj.CalcAdmixedGRMMatrix(pop_af(:,[1,snp_ind]),pop);
causal_grm = obj.CalcGRMMatrix();
geno_d = obj.UnpackGeno();
parfor i=1:num_traits
    random_effect = geno_d(:,2:end)*normrnd(0,1,[num_causal_snps 1]);
    env_effect = zscore(normrnd(0,1,[num_samples,1]));
    gene_effect = zscore(random_effect)*0.7+zscore(geno_d(:,1))*0.3;
    env_effect = env_effect-dot(env_effect,gene_effect)/dot(gene_effect,gene_effect).*gene_effect;
    sim_traits(:,i) = gene_effect+zscore(env_effect)*sqrt(0.42); 
end
[~,h2_r_1] = CalcUniHeritability(sim_traits,[],[],grm,8);
[h2_f_2,h2_r_2] = CalcUniHeritability(sim_traits,[],scores,grm,8);
[h2_f,h2_r] = CalcUniHeritability(sim_traits,[],scores,true_grm,8);
T = table(h2_r_1',h2_f_2',h2_r_2',h2_f',h2_r');
writetable(T,'fam_results.txt');
[~,h2_r_1] = CalcUniHeritability(sim_traits,[],[],causal_grm,8);
[h2_f_2,h2_r_2] = CalcUniHeritability(sim_traits,[],scores,causal_grm,8);
[h2_f,h2_r] = CalcUniHeritability(sim_traits,[],scores,true_causal_grm,8);
T = table(h2_r_1',h2_f_2',h2_r_2',h2_f',h2_r');
writetable(T,'fam_causal_results.txt');