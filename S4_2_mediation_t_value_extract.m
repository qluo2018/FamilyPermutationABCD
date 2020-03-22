%Extract average mediation t value using 4mm ROI centered by AHBA samples' MNI in
%left hemisphere

%% written by Shen Chun, cshen17@fudan.edu.cn
%% reviewed by Dr Qiang Luo, qluo@fudan.edu.cn
%% released on 21 Mar 2020
%% please cite: Shen, et al. Biological Psychiatry 2020

load('ROI_index_80.mat')
%code could be found in S4_1
%load ABCD mask
load('DATA.mat','mask_index')

%whole brain mediation map, a*b t value
ab_1 = spm_vol('path_AB_t.nii');
ab_t = spm_read_vols(ab_1);
ab_tvalue = ab_t(mask_index);

cd('roi_mask')
Index_cortL = Ncort_l_80(:,1);
tvalue_cort_l = zeros(length(Index_cortL),1);
for i = 1:length(Index_cortL)
    ROIfile = strcat('sample_',num2str(Index_cortL(i)),'_sphere.nii');
    v = spm_vol(ROIfile);
    v1 = spm_read_vols(v);
    v2 = find(v1~=0);
    [is,Locb] = ismember(v2,origin_index);%some voxels are not in mask_index
    Locb_use = Locb(is);
    tvalue_sphere =ab_tvalue(Locb_use);
    tvalue_cort_l(i) = nanmean(tvalue_sphere);
end

Index_subL_use = Nsub_l_80_use;
tvalue_sub_l = zeros(length(Index_subL_use),1);
for i = 1:length(Index_subL_use)
    ROIfile = strcat('sample_',num2str(Index_subL_use(i)),'_sphere.nii');
    v = spm_vol(ROIfile);
    v1 = spm_read_vols(v);
    v2 = find(v1~=0);
    [is,Locb] = ismember(v2,origin_index);%some voxels are not in mask_index
    Locb_use = Locb(is);
    tvalue_sphere =ab_tvalue(Locb_use);
    tvalue_sub_l(i) = nanmean(tvalue_sphere);
end

save ab_tvalue_80ROI tvalue_cort_l tvalue_sub_l;
