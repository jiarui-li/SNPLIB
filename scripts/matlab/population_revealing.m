clear all
%% Simulation settings
Fst = [0.001,0.005,0.01,0.05,0.1];
num_generations = 20;
effective_sample_size = floor(num_generations./2./(1-exp(-Fst)));
num_snps = 10000;
num_samples = 1000;
pca_corr = zeros(100,5);
mds_corr = zeros(100,5);
sugibs_corr = zeros(100,5);
for k=1:5    
    tic
    for i=1:100
        %% Simulating Population Structure
        aaf = unifrnd(0.1,0.9,[1,num_snps]);
        pop_af = UpdateAf_(aaf,2,num_generations,effective_sample_size(k))';
        pop = zeros(num_samples, 2);
        pop(:,1) = betarnd(0.5,0.5,[num_samples 1]);
        pop(:,2) = 1- pop(:,1);
        af =  pop*pop_af;
        %% Simulating Genotypes and scores
        obj = SNPLIB();
        obj.nThreads = 4;
        obj.GenerateIndividuals(af);
        pca = obj.CalcPCAScores(1);
        sugibs = obj.CalcSUGIBSScores(1);
        ibs = obj.CalcIBSMatrix();
        mds = cmdscale(ibs,1);
        pca_corr(i,k) = abs(corr(pca,pop(:,1)));
        sugibs_corr(i,k) = abs(corr(sugibs,pop(:,1)));
        mds_corr(i,k) = abs(corr(mds,pop(:,1)));
    end
    toc
end