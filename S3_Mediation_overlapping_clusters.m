%%%%%%%%%%%%%%
%% Step 3: mediation analysis
%%%%%%%%%%%%%%
%% Step 3.1: mediation of the significant clusters
% Using the mediation toolbox developed by Tor Wager,
%http://wagerlab.colorado.edu/tools

%% written by Shen Chun, cshen17@fudan.edu.cn
%% reviewed by Dr Qiang Luo, qluo@fudan.edu.cn
%% released on 21 Mar 2020
%% please cite: Shen, et al. Biological Psychiatry 2020

addpath(genpath('/home1/shenchun/Documents/toolbox/function/')); %3M package

load('DATA.mat'); 
% Descript of DATA.mat:
% GMV --- 3515-by-529551 GMV of each voxel from each subject
% ID --- subject ID
% cov_comb  ---- covariates consideered: using puberty, handedness (2 dummy variables), 
%               household income, age, site (20 dummy variables), TIV, race (3 dummy
%               variables), parental education, BMI, headmotion (provided by
%               ABCD)
% dsm-adhd  ---- CBCL_adhd
% dysomnia  ---- Sleep disturbances rated using the parent Sleep Disturbance
%       Scale for Children (60) were further summarized into two dimensions: 
%       dyssomnias (disorders of initiating and maintaining sleep, sleep 
%       breathing disorders and disorders of excessive somnolence) and 
%       parasomnias (disorders of arousal, sleep-wake transition disorders 
%       and sleep hyperhidrosis) (61). 
% mask_index --- defined as GMV  > 10%

% Extract overlapping clusters
dysomnia_m = spm_vol('dysomnia_per5000_0001_sigclusters_mask.nii');
dysomnia_mask = spm_read_vols(dysomnia_m);

origin_1 = reshape(mask_index,[1,121*145*121]);
origin_index = find(origin_1>0);

dsm_m = spm_vol('dsm_adhd_per5000_0001_sigclusters_mask.nii');
dsm_mask = spm_read_vols(dsm_m);

overlap_dsm_slp = zeros(size(mask_index));
overlap_dsm_slp(dysomnia_mask==1 & dsm_mask ==1)=1;
overlap_dsm_slp_index = find(overlap_dsm_slp>0);

[~,Locb_dsm_slp] = ismember(overlap_dsm_slp_index,origin_index);
GMV_overlap_dsm_slp = GMV(:,Locb_dsm_slp);
GMV_overlap_dsm_slp_mean = mean(GMV_overlap_dsm_slp,2);

% Mediation: brain-->dsm adhd -->sleep
[paths, stats] = mediation(GMV_overlap_dsm_slp_mean,dysomnia,dsm_adhd,...
    'boottop','bootsamples',10000,'covs',cov_comb);
save med_overlap paths stats;
%% Step 3.2: exploratory mediation in the whole brain
% IMPORTANT: we had to run this program by blocks of voxels, becasue the data were too
% large to be processed in a parallel computing setting 
% Load data
load('DATA.mat')

X = GMV;
Y= dysomnia; 
M = dsm_adhd;
cov_all = cov_comb;

% Mediation analysis
index_1 = 1:1000:size(X,2);
index_2 = [1000:1000:size(X,2) size(X,2)];

for j=1:length(index_1)
    j
    X_select = X(:,index_1(j):index_2(j));
    
    paths_1 = cell(size(X_select,2),1);
    stats_1 = cell(size(X_select,2),1);
    ste_1 = cell(size(X_select,2),1);
    parfor i=1:size(X_select,2)
        [paths_2, stats_2] = mediation(X_select(:,i),Y,M,'boottop','bootsamples',3000,'covs',cov_all);
        
        paths_1{i} = stats_2.paths;
        stats_1{i} = stats_2.p;
        ste_1{i} = stats_2.ste;
    end    
    eval(['save result_' num2str(j) '.mat paths_1 stats_1 ste_1'])
end

path_beta_all = [];
path_pval_all = [];
path_ste_all = [];
for i=1:length(index_1)
    eval(['load result_' num2str(i) '.mat'])
    path_beta_all = [path_beta_all; paths_1];
    path_pval_all = [path_pval_all; stats_1];
    path_ste_all = [path_ste_all;ste_1];
end

path_beta_all = cell2mat(path_beta_all);
path_pval_all = cell2mat(path_pval_all);
path_ste_all = cell2mat(path_ste_all);

name_path = {'path_A'; 'path_B'; 'path_C1'; 'path_C'; 'path_AB'};
table_path_pval_all = array2table(path_pval_all, 'VariableNames',name_path);
table_path_beta_all = array2table(path_beta_all, 'VariableNames',name_path);
table_path_ste_all = array2table(path_ste_all, 'VariableNames',name_path);
save wholeBrainMediation_dysomnia table_path_pval_all table_path_beta_all table_path_ste_all;
