%Whole-brain voxel-wise mediation using the mediation toolbox developed by 
%Tor Wager,
%http://wagerlab.colorado.edu/tools

addpath(genpath('/home1/shenchun/Documents/toolbox/Mediation_3M/'));

%load data
load('/share/home1/shenchun/Documents/Projects/sleep/DATA_brain_0426.mat','GMV','dsm_adhd','cov_comb')
load('/share/home1/shenchun/Documents/Projects/sleep/DATA_beh_0426.mat','dysomnia')
load('/share/home1/shenchun/Documents/DATA/ABCD_VBM/smwp1_data_masked.mat','mask_index')

X = GMV;
Y= dysomnia; 
M = dsm_adhd;
cov_all = cov_comb;

%start mediation
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
save path_pval_beta_all_dysomnia table_path_pval_all table_path_beta_all table_path_ste_all table_path_z_all;
