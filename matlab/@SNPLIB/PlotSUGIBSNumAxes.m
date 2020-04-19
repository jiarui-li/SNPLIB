function [D_obj,d_sim_upper] = PlotSUGIBSNumAxes(obj, num_simulations)
%PLOTSUGIBSNUMAXES Summary of this function goes here
%   Detailed explanation goes here
if nargin<2
    num_simulations = 100;
end
ibs = obj.CalcIBSMatrix();
d = sum(ibs);
d = diag(d.^(-0.5));
ugrm = CalcUGRMMatrix_(obj.GENO,obj.nSamples,obj.nThreads);
I = d*ugrm*d;
D = eig(I);
D = sort(D,'descend');
D_obj = D(2:end);
D_obj = D_obj - mean(D_obj);
af = obj.CalcAlleleFrequency();
D_sim = zeros(obj.nSamples-1,num_simulations);
n = obj.nSamples;
parfor i=1:num_simulations
    GENO = GenerateSimGENO_(af, n);
    ibs = CalcIBSMatrix_(GENO,n,1);
    ugrm = CalcUGRMMatrix_(GENO,n,1);
    d = sum(ibs);
    d = diag(d.^(-0.5));
    I = d*ugrm*d;
    D = eig(I);
    D = sort(D,'descend');
    D_sim(:,i) = D(2:end);
end
D_sim = bsxfun(@minus,D_sim,mean(D_sim));
d_sim = median(D_sim,2);
b = robustfit(d_sim,D_obj,'bisquare',4.685,'off');
D_sim = b.*D_sim;
d_sim_upper = quantile(D_sim,0.95,2);
plot([D_obj(1:30),d_sim_upper(1:30)],'.')
end

