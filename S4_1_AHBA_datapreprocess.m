%AHBA data preprocessing, code adapted from Arnatkeviciute et al.(2019)
%Reference: Arnatkevic Iute A, Fulcher BD, Fornito A (2019): 
%A practical guide to linking brain-wide gene expression and 
%neuroimaging data. Neuroimage 189:353-367.

%% written by Shen Chun, cshen17@fudan.edu.cn
%% reviewed by Dr Qiang Luo, qluo@fudan.edu.cn
%% released on 21 Mar 2020
%% please cite: Shen, et al. Biological Psychiatry 2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Step 1: Reannotation: Load AHBA raw data and reannotated data (Method S3, processed by Xingzhong Zhao)
% The reannotation was done by using the tool provided by 
% https://github.com/BMHLab/AHBAprocessing/blob/master/code/dataProcessing/README_Reannotator.txt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

remap = readtable('AHBRemappingFile.txt');
probeID_n = table2array(remap(:,1));
probeName_n = table2array(remap(:,2));
geneName_n = table2array(remap(:,3));
entrezID_n = str2double(table2array(remap(:,4)));

%Load allen brain data
cd('AHBA_data')
% https://human.brain-map.org
% 6 donors, each donor has a data folder

%AHBA proble info
FileProbes = './normalized_microarray_donor9861/Probes.xlsx';
ProbeTable = readtable(FileProbes);
ProbeID = ProbeTable.probe_id;
ProbeName =  ProbeTable.probe_name;
GeneSymbol = ProbeTable.gene_symbol;
EntrezID = ProbeTable.entrez_id;

[~,ia,ib] = intersect(ProbeName, probeName_n);%46051,need to match by probe ID
EntrezID = entrezID_n(ib);
ProbeName = ProbeName(ia);
GeneSymbol = geneName_n(ib);

headerdata = {'Expression', 'MNIcoordinates', 'StructureName', 'MRIvoxCoordinates', 'Noise', 'SampleID', 'WellID'};
headerprobe = { 'ProbeName', 'EntrezID', 'GeneSymbol'};

%Go to each subject's directory and take the data
donor = [9861,10021,12876,14380,15496,15697];

for subj=1:6
    folder = sprintf('normalized_microarray_donor%d', donor(subj));
    cd(folder);
    %%load information specific for each subject
    FileMicroarray = 'MicroarrayExpression.csv';
    FileAnnot = 'SampleAnnot.xlsx';
    FileNoise = 'PACall.csv';
    Expression = csvread(FileMicroarray);
    noise = csvread(FileNoise);
    noise = noise(ia,:);

    Expression(:,1) = [];                         % exclude probe IDs from expression matrix
    Expression = Expression(ia,:);
    [~,~,SlabType] = xlsread(FileAnnot, 'D:D');
    SlabType(any(cellfun(@(x) any(isnan(x)), SlabType),2),:) = [];

    [~,~, StructureName] = xlsread(FileAnnot, 'F:F');
    StructureName(any(cellfun(@(x) any(isnan(x)), StructureName),2),:) = [];

    SlabType(1) = [];                           % remove headline
    StructureName(1) = [];                      % remove headline
    MMcoordinates = xlsread(FileAnnot, 'K:M');
    MNI_index = round(MMcoordinates,0); %round to integer
    MRIvoxCoordinates = xlsread(FileAnnot, 'H:J');
    SampleID = xlsread(FileAnnot, 'A:A');
    WellID = xlsread(FileAnnot, 'C:C');
    [~,probeList] = intersect(noise(:,1),ProbeID, 'stable');%make sure probe ID sequence is the same
    noise = noise(probeList,2:end);

    % for nan columns
    % keep only existing expression values
    Expression = Expression(:,all(~isnan(Expression)));
    % keep only existing coordinates
    MNI_index = MNI_index(all(~isnan(MNI_index),2),:); % for nan rows
    MRIvoxCoordinates = MRIvoxCoordinates(all(~isnan(MRIvoxCoordinates),2),:); % for nan rows
    SampleID = SampleID(all(~isnan(SampleID),2),:); % for nan rows
    WellID = WellID(all(~isnan(WellID),2),:); % for nan rows
    noise = noise(:,all(~isnan(noise)));
    % keep only existing structure names
    StructureName(cellfun(@(StructureName) any(isnan(StructureName)),StructureName)) = [];

    % assign output to Data cell;
    Data{subj,1} = Expression;
    Data{subj,2} = MNI_index;
    Data{subj,3} = StructureName;
    Data{subj,4} = MRIvoxCoordinates;
    Data{subj,5} = noise;
    Data{subj,6} = SampleID;
    Data{subj,7} = WellID;
    cd('..')
end

DataTable = dataset({Data, headerdata{:}});

%remove probes with missing entrezIDs
for s=1:6
    DataTable.Expression{s,1}(isnan(EntrezID),:) = [];
    DataTable.Noise{s,1}(isnan(EntrezID),:) = [];
end

ProbeName(isnan(EntrezID)) = [];
GeneSymbol(isnan(EntrezID)) = [];
EntrezID(isnan(EntrezID)) = [];

fprintf(1,'%d unique genes\n', length(unique(EntrezID)))

% Assign ProbeIDs, EntrezIDs and ProbeNames to Data cell.
DataProbe{1,1} = ProbeName;
DataProbe{1,2} = EntrezID;
DataProbe{1,3} = GeneSymbol;
DataTableProbe = dataset({DataProbe, headerprobe{:}});

% Combine data for all subjects
Expressionall = horzcat(DataTable{1,1}, DataTable{2,1}, DataTable{3,1}, DataTable{4,1}, DataTable{5,1}, DataTable{6,1});
Coordinatesall = vertcat(DataTable{1,2}, DataTable{2,2}, DataTable{3,2}, DataTable{4,2}, DataTable{5,2}, DataTable{6,2});
StructureNamesall = vertcat(DataTable{1,3}, DataTable{2,3}, DataTable{3,3}, DataTable{4,3}, DataTable{5,3}, DataTable{6,3});
MRIvoxCoordinatesAll = vertcat(DataTable{1,4}, DataTable{2,4}, DataTable{3,4}, DataTable{4,4}, DataTable{5,4}, DataTable{6,4});
noiseall = horzcat(DataTable{1,5}, DataTable{2,5}, DataTable{3,5}, DataTable{4,5}, DataTable{5,5}, DataTable{6,5});

save AHBA_data_reanote DataTable DataTableProbe Expressionall Coordinatesall StructureNamesall MRIvoxCoordinatesAll noiseall -v7.3;

%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 2. Data filtering: exclude probes that don't exceed the background in 
% at least 50% of all cortical and subcortical samples across all subjects
%%%%%%%%%%%%%%%%%%%%%%%%

signalThreshold = 0.5;

signalLevel = sum(noiseall,2)./size(noiseall,2);%1: above background
indKeepProbes = find(signalLevel>=signalThreshold);

% remove selected probes from data and perform the following analysis using
% the non-noisy probes only
% original sequence
ProbeName = DataTableProbe.ProbeName{1,1}(indKeepProbes);
ProbeID = DataTableProbe.ProbeID{1,1}(indKeepProbes);
EntrezID = DataTableProbe.EntrezID{1,1}(indKeepProbes);
GeneSymbol = DataTableProbe.GeneSymbol{1,1}(indKeepProbes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 3: Probe selection: Using the mean expressions among the probes for the same gene
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sort by EntrezID
[~, sortALL] = sort(EntrezID);
EntrezID_s = EntrezID(sortALL);
GeneSymbol_s =  GeneSymbol(sortALL);

nSub = 6;
Uniq = unique(EntrezID_s);
expressionMean = cell(6,1);
for subj = 1:nSub
    expression = (DataTable.Expression{subj}(indKeepProbes,:))';
    for k=1:length(Uniq)
        % find indexes for repeating entrexIDs
        indRepEntrezIDs = find(EntrezID==(Uniq(k)));
        expRepEntrezIDs = expression(:,indRepEntrezIDs);
        if length(indRepEntrezIDs) >=2  
            expressionMean{subj}(:,k) = mean(expRepEntrezIDs,2);%column sorted by EntrezID_s
        else
            expressionMean{subj}(:,k) = expRepEntrezIDs;
        end   
    end
end

[c,ia,ib] = intersect(Uniq,EntrezID_s);
probeInformation.EntrezID = Uniq; 
probeInformation.GeneSymbol = GeneSymbol_s(ib);

sampleInfo = cell(6,1);
for subject=1:6
    % combine sample information variables to a structure.
    SampleInformation.StructureNames = DataTable.StructureName{subject,1};
    SampleInformation.StructureID = DataTable.SampleID{subject,1};
    SampleInformation.MNICoordinates = DataTable.MNIcoordinates{subject,1};
    SampleInformation.MRIvoxCoordinates = DataTable.MRIvoxCoordinates{subject,1};
    sampleInfo{subject} = SampleInformation;
end

save AHBA_Mean_reanote expressionMean probeInformation sampleInfo;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 4. Separating the samples into cortical and subcortical 
% Define the cortical and the subcortical ROIs based on the MNI coordinates
% of the AHBA samples.  
% When more than 80% of a small ball centered at the location of each sample 
% (i.e. 4mm sphere ROI) overlapped with the Harvard mask and the ABCD mask,
% this sample was classified as either a cortical or a subcortical sample.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load HarvardOxford cortical atlas
v1 = spm_vol('Reslice_HarvardOxford-cort.nii');
cort = spm_read_vols(v1);
cort_index = find(cort~=0);
cort_l = cort;%left hemisphere
cort_l(mod(cort_l,2)==1) = 0;%left is even, odd is right
cort_l_index = find(cort_l~=0);

%Load HarvardOxford subcortical atlas
v2 = spm_vol('Reslice_HarvardOxford-sub.nii');
sub = spm_read_vols(v2);
sub_index = find(sub~=0);
sub_l = sub;%left hemisphere
sub_l(sub_l ==1116|sub_l ==2049|sub_l ==3050|sub_l ==4051|sub_l ==5052|sub_l ==6053|sub_l ==7054|sub_l ==8058) = 0;
sub_l_index = find(sub_l~=0);

%ABCD GMV mask
load('DATA.mat','mask_index')
origin_1 = reshape(mask_index,[1,121*145*121]);
origin_index = find(origin_1>0);
origin_index = origin_index';

% to have the header information of the ABCD sMRI images
ab_1 = spm_vol('dysomnia_per5000_0001_sigclusters_mask.nii'); 

%MNI coordinates of AHBA samples
load('AHBA_data_reanote.mat','Coordinatesall')

len = size(Coordinatesall,1);
MNI_x = Coordinatesall(:,1);
MNI_y = Coordinatesall(:,2);
MNI_z = Coordinatesall(:,3);

%create ROI mask using FSL
%make a new folder to save ROI masks
mkdir('roi_mask')
tvalue_ROI = [];
for i = 1:len
    %transfer MNI coordinates to voxel coordinates
    voxel_index = inv(ab_1.mat)*[MNI_x(i);MNI_y(i);MNI_z(i);1];
    x = round(voxel_index(1),0);
    y = round(voxel_index(2),0);
    z = round(voxel_index(3),0);
    %input1 = '/share/home1/shenchun/Documents/Projects/sleep/FINAL_0426/3_mediation/3_Wholebrain_Mediation/path_AB_t.nii';
    output1 = strcat('sample_',num2str(i),'_point');
    %create single point mask
    unix(['fslmaths mask_orig.nii -mul 0 -add 1 -roi ',...
        num2str(x),' 1 ',num2str(y),' 1 ',num2str(z),' 1 0 1 ',output1,' -odt float']);
    %create sphere
    output2 = strcat('sample_',num2str(i),'_sphere');
    unix(['fslmaths ',output1,' -kernel sphere 4 -fmean ',output2,' -odt float']);
    %gunzip
    gunzip([output2,'.nii.gz']);
    spherefile = strcat('sample_',num2str(i),'_sphere.nii');
end

%ROI in left hemisphere
Ncort_l = zeros(len,1);%left
Nsub_l = zeros(len,1);

cd('roi_mask')
for i = 1:len
    sample = strcat('sample_',num2str(i),'_sphere.nii');
    v3 = spm_vol(sample);
    vv3 = spm_read_vols(v3);
    index = find(vv3~=0);
    [is1,~] = ismember(index,cort_l_index);
    [is2,~] = ismember(index,sub_l_index);
    Ncort_l(i) = sum(is1,1)./length(index);
    Nsub_l(i) = sum(is2,1)./length(index);
end

%80% overlap with Harvard mask
Ncort_l_80 = [find(Ncort_l>=0.8),Ncort_l(Ncort_l>=0.8)];%n=784
Nsub_l_80 = [find(Nsub_l>=0.8),Nsub_l(Nsub_l>=0.8)];%n=437

%80% overlap with ABCD GMV mask
%cortical ROIs,left
Index_cortL = Ncort_l_80(:,1);
cort_per = zeros(length(Index_cortL),1);
for i = 1:length(Index_cortL)
    ROIfile = strcat('sample_',num2str(Index_cortL(i)),'_sphere.nii');
    v = spm_vol(ROIfile);
    v1 = spm_read_vols(v);
    v2 = find(v1~=0);
    [is,~] = ismember(v2,origin_index);%some voxels are not in mask_index
    cort_per(i) = sum(is,1)./length(v2);%most of ROIs in mask_index
end

%subcortical ROIs,left
Index_subL = Nsub_l_80(:,1);
sub_per = zeros(length(Index_subL),1);
for i = 1:length(Index_subL)
    ROIfile = strcat('sample_',num2str(Index_subL(i)),'_sphere.nii');
    v = spm_vol(ROIfile);
    v1 = spm_read_vols(v);
    v2 = find(v1~=0);
    [is,~] = ismember(v2,origin_index);%some voxels are not in mask_index
    sub_per(i) = sum(is,1)./length(v2);%some ROIs not in mask index
end
Nsub_l_80_use = Nsub_l_80(find(sub_per>=0.8),1);%n=182
save ROI_index_80 Ncort_l_80 Nsub_l_80_use;
%%
%extract cortical, subcortical gene expression data, 
%based on HarOx atlas, 80% ~NaN, left hemisphere
load('AHBA_Mean_reanote.mat')
load('ROI_index_80.mat')
%%
cort_expMeanRaw = cell(6,1);
sub_expMeanRaw = cell(6,1);

Index_cortL = Ncort_l_80(:,1);
Index_subL_use = Nsub_l_80_use;

%samples by individual
%cortical
for i = 1:length(Index_cortL)
    if Index_cortL(i)<947
        Index_cortL(i,2) = 1;
    elseif Index_cortL(i)>946 && Index_cortL(i)<1840
        Index_cortL(i,2) = 2;
    elseif Index_cortL(i)>1839 && Index_cortL(i)<2203
        Index_cortL(i,2) = 3;
    elseif Index_cortL(i)>2202 && Index_cortL(i)<2732
        Index_cortL(i,2) = 4;
    elseif Index_cortL(i)>2731 && Index_cortL(i)<3202
        Index_cortL(i,2) = 5;
    elseif Index_cortL(i)>3201
        Index_cortL(i,2) = 6;
    end
end
tabulate(Index_cortL(:,2))

%subcortical
for i = 1:length(Index_subL_use)
    if Index_subL_use(i)<947
        Index_subL_use(i,2) = 1;
    elseif Index_subL_use(i)>946 && Index_subL_use(i)<1840
        Index_subL_use(i,2) = 2;
    elseif Index_subL_use(i)>1839 && Index_subL_use(i)<2203
        Index_subL_use(i,2) = 3;
    elseif Index_subL_use(i)>2202 && Index_subL_use(i)<2732
        Index_subL_use(i,2) = 4;
    elseif Index_subL_use(i)>2731 && Index_subL_use(i)<3202
        Index_subL_use(i,2) = 5;
    elseif Index_subL_use(i)>3201
        Index_subL_use(i,2) = 6;
    end
end
tabulate(Index_subL_use(:,2))

%extract cortical and subcortical raw data
%cortical
for subject=1:6
     if subject ==1
         exp_raw = expressionMean{subject}(Index_cortL(Index_cortL(:,2)==subject,1),:); 
     elseif subject ==2
         index = Index_cortL(Index_cortL(:,2)==subject,1)-946;
         exp_raw = expressionMean{subject}(index,:);
     elseif subject ==3
         index = Index_cortL(Index_cortL(:,2)==subject,1)-1839;
         exp_raw = expressionMean{subject}(index,:);
     elseif subject ==4
         index = Index_cortL(Index_cortL(:,2)==subject,1)-2202;
         exp_raw = expressionMean{subject}(index,:);
     elseif subject ==5
         index = Index_cortL(Index_cortL(:,2)==subject,1)-2731;
         exp_raw = expressionMean{subject}(index,:);
     elseif subject ==6
         index = Index_cortL(Index_cortL(:,2)==subject,1)-3201;
         exp_raw = expressionMean{subject}(index,:);
     end
     cort_expMeanRaw{subject} = exp_raw;
end

%subcortical
exp_raw = [];
for subject=1:6
     if subject ==1
         exp_raw = expressionMean{subject}(Index_subL_use(Index_subL_use(:,2)==subject,1),:); 
     elseif subject ==2
         index = Index_subL_use(Index_subL_use(:,2)==subject,1)-946;
         exp_raw = expressionMean{subject}(index,:);
     elseif subject ==3
         index = Index_subL_use(Index_subL_use(:,2)==subject,1)-1839;
         exp_raw = expressionMean{subject}(index,:);
     elseif subject ==4
         index = Index_subL_use(Index_subL_use(:,2)==subject,1)-2202;
         exp_raw = expressionMean{subject}(index,:);
     elseif subject ==5
         index = Index_subL_use(Index_subL_use(:,2)==subject,1)-2731;
         exp_raw = expressionMean{subject}(index,:);
     elseif subject ==6
         index = Index_subL_use(Index_subL_use(:,2)==subject,1)-3201;
         exp_raw = expressionMean{subject}(index,:);
     end
     sub_expMeanRaw{subject} = exp_raw;
end 
save exp_Mean_80ROI_raw_reanote Index_cortL Index_subL_use cort_expMeanRaw sub_expMeanRaw;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 5: normalization within donors
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BF_NormalizeMatrix function written by Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
addpath(genpath('/home1/shenchun/Documents/toolbox/function'))

%scaledRobustSigmoid
cort_expMeanScaled = cell(6,1);
sub_expMeanScaled = cell(6,1);
cort_badNorm = cell(6,1);%badly-normalised genes
sub_badNorm = cell(6,1);

%first normalized with sample across gene, then normalized across sample
for subject = 1:6
    exp_raw = cort_expMeanRaw{subject};
    exp_norm1 = BF_NormalizeMatrix(exp_raw','scaledRobustSigmoid')';
    exp_norm2 = BF_NormalizeMatrix(exp_norm1,'scaledRobustSigmoid');
    temp  = find(any(isnan(exp_norm2),1));
    cort_expMeanScaled{subject} = exp_norm2;
    cort_badNorm{subject} = temp;
end

for subject = 1:6
    exp_raw_s = sub_expMeanRaw{subject};
    exp_norm_s1 = BF_NormalizeMatrix(exp_raw_s','scaledRobustSigmoid')';
    exp_norm_s2 = BF_NormalizeMatrix(exp_norm_s1,'scaledRobustSigmoid');
    temp2  = find(any(isnan(exp_norm_s2),1));
    sub_expMeanScaled{subject} = exp_norm_s2;
    sub_badNorm{subject} = temp2;
end

badNorm = horzcat(cort_badNorm{1},cort_badNorm{2},cort_badNorm{3},...
    cort_badNorm{4},cort_badNorm{5},cort_badNorm{6},...
    sub_badNorm{1},sub_badNorm{2},sub_badNorm{3},...
    sub_badNorm{4},sub_badNorm{5},sub_badNorm{6});
ubadNorm = unique(badNorm);

cort_l_expMS = vertcat(cort_expMeanScaled{1},cort_expMeanScaled{2},cort_expMeanScaled{3},...
    cort_expMeanScaled{4},cort_expMeanScaled{5},cort_expMeanScaled{6});
sub_l_expMS = vertcat(sub_expMeanScaled{1},sub_expMeanScaled{2},sub_expMeanScaled{3},...
    sub_expMeanScaled{4},sub_expMeanScaled{5},sub_expMeanScaled{6});

%remove the bad normalized genes
cort_l_expMS(:,ubadNorm) = [];
sub_l_expMS(:,ubadNorm) = [];

load('AHBA_Mean_reanote.mat','probeInformation')
genes=probeInformation.EntrezID; 
geneSymbol = probeInformation.GeneSymbol;

genes(ubadNorm) = [];
geneSymbol(ubadNorm) = [];

save AHBA_Mean_scaled_reanote cort_expMeanScaled cort_l_expMS...
    genes geneSymbol sub_expMeanScaled sub_l_expMS;