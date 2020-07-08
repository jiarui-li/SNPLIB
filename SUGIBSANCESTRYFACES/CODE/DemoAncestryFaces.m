%% SCRIPT TO DEMONSTRATE THE RECONSTRUCTION OF ANCESTRY FACES
close all;clear all;
restoredefaultpath;
localpath = '.../SUGIBSANCESTRYFACES/';% provide your local path to the code and data
datapath = [localpath 'DATA/'];% define data path
addpath(genpath([localpath 'CODE/']));% add code to matlab environment
%% LOADING DATA
load([datapath 'AVGscan.mat']);% Loading the average ancestry face as reference face, can be viewed using viewer(AVGscan);
load([datapath 'AncientInfo.mat']);% Loading information on the ancient DNA profiles
load([datapath 'ONEKG.mat']);% Loading information on the 1000 genome project
load([datapath 'PLSR.mat']);% Loading the PLSR (BETA) mordel extracted from the PSU cohort as projected in 1KG ancestry space;
%% SETTING OTHER COVARIATES
age = 30;% years
sex = 1.5;% 1 = males, 2 = females
height = 165;% cm
weight = 75;%kg
BASEX = [age sex height weight zeros(1,8)];
%% RECONSTRUCTING ANCESTRY FACES FOR THE ANCIENT DNA PROFILES-> FIGURE 6
n = length(ancient.labels);
for i=1:1:n
    % i=1;
    % reconstruct face
    X = BASEX;
    X(2) = ancient.sex(i);% change to approriate sex value;
    X(5:end) = ancient.ancestry(i,:);% change to approriate ancestry scores
    Y = [1 X]*BETA;% apply PLSR model
    Y = reshape(Y(:),3,AVGscan.nVertices);% get result into shape vertices
    scan = clone(AVGscan);
    scan.Vertices = Y';
    v = viewer(scan);
    v.SceneLightVisible = true;
    v.SceneLightLinked = true;
    title(scan.RenderAxes,ancient.labels{i},'Color',[1 1 1]);
end
%% RECONSTRUCTING ANCESTRY FACES TO ILLUSTRATE GENOMIC VARIATION -> FIGURE 5
% select axis to illustrate
pcax = 1;% first axis selected
refpc = ONEKG.ancestry(:,pcax);meanpc = mean(refpc);stdpc = std(refpc);% collecting reference information (mean and standard deveiation) based on the 1KG project data
% reconstruct face along positive side of the axis
X = BASEX;
X(4+pcax) = meanpc+3*stdpc;% change to approriate ancestry scores
Y = [1 X]*BETA;% apply PLSR model
Y = reshape(Y(:),3,AVGscan.nVertices);% get result into shape vertices
scan = clone(AVGscan);
scan.Vertices = Y';
v = viewer(scan);
v.SceneLightVisible = true;
v.SceneLightLinked = true;
title(scan.RenderAxes,'Plus 3std','Color',[1 1 1]);
% reconstruct face along positive side of the axis
X = BASEX;
X(4+pcax) = meanpc-3*stdpc;% change to approriate ancestry scores
Y = [1 X]*BETA;% apply PLSR model
Y = reshape(Y(:),3,AVGscan.nVertices);% get result into shape vertices
scan = clone(AVGscan);
scan.Vertices = Y';
v = viewer(scan);
v.SceneLightVisible = true;
v.SceneLightLinked = true;
title(scan.RenderAxes,'Minus 3std','Color',[1 1 1]);
%%