%% Example for running SigProfilerSingleSample 
%% The example examines 10 biliary adenocarcinoma whole-genomes from ICGC PCAWG and
%% assigns mutational signatures using upcoming release of PCAWG mutational signatures
%% Clearing all data
close all;
clearvars;
clc;
addpath('source/');

%% Starting the default cluster
if ( ~isempty(gcp('nocreate')) )
    delete(gcp);
end
c = parcluster;
job = parpool(c);

%% Analysis of all signatures in individual samples
tic
    analysis_individual_samples('input/Consensus_subs_mutational_signatures.mat', ... % set of signatures
                                'input/signatures_in_samples_and_cancer_types.mat', ... % set of signatures in samples
                                'input/Biliary-AdenoCA_example_cancer_samples.mat', ... % set of individual samples for examination
                                'output/', ... % output folder
                                'signatures_in_Biliary-AdenoCA_example_cancer_samples.mat', ... % output file
                                1, ... % boolean variable indicating whether to use rules or not (1==use rules; 0==do not use rules)
                                [1 5], ... % IDs of signatures to be included in all samples regardless of rules or sparsity  
                                {[2 17], [7 8 9 10], [13 14], [21 22]}); % connected signatures (e.g., signatures SBS-2 and SBS-13)
toc