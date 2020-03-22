%% To define the exchangable blocks by using the family relatedness 
% extracted from a questionnaire: acspsw02 in the ABCD Data Release v1.0
%% Input:
% 'familyinfo.mat' ---- 3 variables: ID (N-by-1 cell), famRela(N-by-1 double), famID (N-by-1, double);
% 'ID_3076.mat'    ---- 1 variable: M-by-1 cell.
%% Output: 
% The multi-level blocks was saved in 'EB_abcd_3076.mat' with two variables:
%                    B --- as an input to palm_quickperms (provided by PALM -- Permutation Analysis of Linear Models); 
%                    ID_3076 --- subject ID. 
% The multi-level block permuted sample in  'crosslag_permorder_3076.mat', namely the 'Pset'.
%% written by Shen Chun, cshen17@fudan.edu.cn
%% reviewed by Dr Qiang Luo, qluo@fudan.edu.cn
%% released on 21 Mar 2020
%% please cite: Shen, et al. Biological Psychiatry 2020

%% NOTE: This code can be adapted to generate the permutation
% sample to perform the multi-level block permutation for the
% participants you have selected from the ABCD cohort.

% Two things need to be dealt with:
% 1. Two data files are neccessary:
% 1.1 'familyinfo.mat' has 3 variables: ID (N-by-1 cell), famRela(N-by-1
% double), famID (N-by-1, double). This information was extrated from the
% acspsw02 questionnaire. The famRela can be 0--singleton, 1--sibling,
% 2--twin, 3--triplet
% 1.2  'ID_3076.mat' has 1 variable: M-by-1 cell. These participant IDs are
% selected according to your criterion
% 2. complex family type
% In the participants selected in our study (Shen and Luo, et al.
% Biological Psychiatry 2020), we had one special family type that had two
% twins and two of three triplets (the third triplet was not seleted into our
% study). Therefore, in this family, the twins were not exchangable with the
% triplets (see the codes in lines 115-116 of this file). 

% Once you have prepared those two data files in point 1, then you can
% check if you have some speical family types as descrbed in point 2. If
% so, you need to deal with them according to the exchangability. 


%kinship based on questionnaire for 3076 crosslag sample
load('familyinfo.mat'); % 3 variables: ID (N-by-1 cell), famRela(N-by-1 double), famID (N-by-1, double); 
                        % N=3515 participants who had the complete behavioural data and qualified VBM data (CAT12, QC > 65) at the baseline
load('ID_3076.mat');    % 3076 of the 3515 participants who had also the complete behaivoural data at the follow-up

[cc,ia,ib] = intersect(ID_3076,ID);
famRela = famRela(ib);
famID = famID(ib);
%tabulate(famRela)
%Value    Count   Percent
%Value    Count   Percent
%      0     2336     75.94%    %singleton
%      1      170      5.53%    %sibling
%      2      564     18.34%    %twin
%      3        6      0.20%    %triplet

F = unique(famID); % 2707 unique families
tt = tabulate(famID);
tt1 = tabulate(tt(:,2));   % the number of participants from each family in this subset of the ABCD cohort
%Value    Count   Percent
%      1     2346     86.67%    
%      2      354     13.08%
%      3        6      0.22%
%      4        1      0.00%

%10 single
famid_s = tt(tt(:,2)==1,1);
idx_s = setdiff(famid_s,famID(famRela==0));
[~,~,ib1] = intersect(idx_s,famID);
%check
%famRela(ib1) %4 sib and 6 twin in 3515 data
famRela(ib1) = 0;
%tabulate(famRela)
%Value    Count   Percent
%      0     2346     76.27%
%      1      166      5.40%
%      2      558     18.14%
%      3        6      0.20%
save kinship_Q_3076 ID_3076 famID famRela;
%% generating the blocks
% Code adapated from hcp2blocks.m by Anderson M. Winkler
% http://brainder.org
% Reference:
% Winkler AM, Webster MA, Vidaurre D, Nichols TE, Smith SM.
% Multi-level block permutation. Neuroimage. 2015;123:253-68.
clear,clc

load('kinship_Q_3076.mat')

N = size(ID_3076,1);

%sibtype
sibtype = zeros(N,1);
for n = 1:N,
    if famRela(n)<2,%single,sibling(not twin)=10
        sibtype(n) = 10;
    elseif famRela(n)==2,
        sibtype(n) = 100;%twin
    elseif famRela(n)==3,
        sibtype(n) = 1000;%triple
    end
end

%famtype
F = unique(famID);
famtype = zeros(N,1);
for f = 1:numel(F),
    fidx = (F(f) == famID);% find subs in same family
    famtype(fidx) = sum(sibtype(fidx));
end
tt1 = tabulate(famtype);
tt1 = tt1(tt1(:,2)~=0,:);
% 6 famtypes [10,20,200,210,400,3000]

%unique subs, two pair of twins in one family
%NDAR_INVEC2CN745 and NDAR_INVV9BM1RB5 are twins based on age
sibtype(strcmp(ID_3076,'NDAR_INVEC2CN745'))=101;
sibtype(strcmp(ID_3076,'NDAR_INVV9BM1RB5'))=101;
%tt2 = tabulate(sibtype);
%tt2 = tt2(tt2(:,2)~=0,:);

[~,idx] = sortrows([famID sibtype]);% ascending sort by famid(1) sibtype(2)
[~,idxback] = sort(idx);
sibtype = sibtype(idx);
famID = famID(idx);
famtype = famtype(idx);

% Now make the blocks for each family
B = cell(numel(F),1);
for f = 1:numel(F),
    fidx = (F(f) == famID);
    ft = famtype(find(fidx,1));% return the first indice =fidx
    if any(ft == [10,20,200,400,3000]),% within-family shuffle
        B{f} = horzcat(famID(fidx),sibtype(fidx));
    else
        B{f} = horzcat(-famID(fidx),sibtype(fidx));
    end
end

% Concatenate all. Prepending the famtype ensures that the
% families of the same type can be shuffled whole-block. Also,
% add column with -1, for within-block at the outermost level
B = horzcat(-ones(N,1),famtype,cell2mat(B));

B = B(idxback,:);
            
save EB_abcd_3076 B ID_3076;

%% calling palm_quickperms to generate the permuted sample
%create permutation order
EBfile = 'EB_abcd_3076.mat';
PSfile = 'crosslag_permorder_3076.mat';
PermSamplePALM(EBfile, PSfile)