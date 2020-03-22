%% family relatedness given by combining the family info and the ABCD genetics
% to have the MZ, DZ, FS, and HS information that was not provided by the questionnaire
% acspsw02.txt in the ABCD study
%% writtent by Dr. Qiang Luo
%% Email: qluo@fudan.edu.cn
%% Release Date: 21 Mar 2020
%% please cite: Shen et al. Biological Psychiatry 2020

%% how to adapt the code to your study:
% First, we loaded the family questionnaire from ABCD acspsw02.txt release
% v1.0 for 4,521 subjects
% Second, we loaded the kinship estimated by KING and the familyinfo.mat, 
% which contains the 3,515 subjects used in our study only. So, we crossed
% these two files by subject ID to generate a kinshipwithin variable with
% kinship information among the 3,515 subjects selected.
% Third, we compared the kinship-estimated family relatedness with the 
% self-reported family relatedness using the cutoffs recommended by the
% KING. We counted the mismatches in each categories:
% >0.35 (MZ), >0.177 (DZ or FS depending on the age difference), >0.0884 (HS).
% Forth, we looked into the mismatches in details. 
% Therefore, to adapt this code for your own study, the 3rd and the 4th 
% steps may need some modifications. 

%% input 
% acspsw02.txt: ---- family qeustionnaire from ABCD release v1.0;  (data for 4,521 subjects)
% kinship_0.0884.txt --- kinship matrix with kinship greater than 0.0884
%                        given by the KING package;
% familyinfo.mat --- family inforamtion extracted from acspsw02.txt: ID, famID, 
%                    famRela(0--singleton, 1--sibling, 2--twin, 3--triplet)
%                    (data for 3,515 with both behaviour and neuroimaging data avaialbe at the baseline)
%          NOTE: 4 famID overlap between single and twin, [584,
%          596,2205,4454], we had re-set those 4 singleton [famID==?? & famRela_final==0] to sibling (1).
% ABCD_release_2.0.1.csv  --- subject ID with gene data releaseed by ABCD v2.0.1  (10,627 subjects)                      
%% output
% imagingInfo2.mat ---- c for subject ID; 
%                       famID_1 for family ID; 
%                       famRela for family relatioinship from the questionannire; 
%                       sibtype  ---- sibling type; 
%                       famtype --- family type

%% Description of the method used here:
% Here, by using the genetic kinship estimated by the KING
% package, we can have more detailed inforamtion about the DZ, MZ, FS and
% HS. However, the kinship was not exactly the same as the family
% relatedness specified in the questionnaire. Therefore, we had done some
% corrections by combing both questionnair information and the kinship. 
% 1) if two subjects had a kinship greater than 0.354, they might
% be MZ. However, if they did not reported the same age, then we had to
% randomly delete one of them from the data set and set the sibling type of
% the remaining subject as 'NS'.
% 2) if two subjeects had a kinship greater than 0.177, they might be DZ
% when they reported the same age, or FS when they reported different ages.
% However, if they reported two different ages with less than 7 months 
% difference, we had to delete randomly one of them from the data set.  
% 3) if two subjects had a kinship greater than 0.0884, they might be HS. 
% 4) if one of twins or siblings was missing the genetic data, we had to
% delete this subject and set the sibling type of the remaining twin or 
% sibling(s) as 'NS'.


%% Summary of the changes made to the sibling type and family ID
% Here, we had 3470 subjects in this file. 
% 45 subject were delted, and among them 42 were deleted owing to the 
% missing genetic data while 3 were deleted owing to different reasons:
% 1) D3H6 had a high kinship with LPC2 (kinship=0.4988), however, these two
% subjects had different ages and could not be a pair of MZ twins.
% 2) 3U3X had a high kinship with 3BZ3 (kinship=0.4988), however, these two
% subjects had different ages and could not be a pair of MZ twins.
% 3) 4K38 had a high kinship with M3Y9 (kinship=0.2474), however, these two
% subjects had different ages and could not be a pair of DZ twins, while
% these two subjects had a 6-month age difference which could not be a
% full-sibling neither.
% The full list was given below:
%     {'NDAR_INVGT92D3H6'}
%     {'NDAR_INVTHEX3U3X'}
%     {'NDAR_INVNNA7JH41'}
%     {'NDAR_INV31JMTLE9'}
%     {'NDAR_INV6NRB61V1'}
%     {'NDAR_INV7CN47GEF'}
%     {'NDAR_INVAEMLZKT4'}
%     {'NDAR_INVAUR9UZH4'}
%     {'NDAR_INVB3438R23'}
%     {'NDAR_INVC6MPJ7ZF'}
%     {'NDAR_INVDLMML76V'}
%     {'NDAR_INVE7CFEXUC'}
%     {'NDAR_INVF19VWPA1'}
%     {'NDAR_INVLDH3YU2R'}
%     {'NDAR_INVMWC4BV86'}
%     {'NDAR_INVVF539UBP'}
%     {'NDAR_INVWDCENDZ9'}
%     {'NDAR_INVWMM5T2PL'}
%     {'NDAR_INVYDEWFLED'}
%     {'NDAR_INVZFYGCHKB'}
%     {'NDAR_INVV7D70A79'}
%     {'NDAR_INVWCV58UKD'}
%     {'NDAR_INVEKN9FWNN'}
%     {'NDAR_INVA6FR087J'}
%     {'NDAR_INVCTH2AK0B'}
%     {'NDAR_INVTWU5RFVX'}
%     {'NDAR_INVH0DGUW54'}
%     {'NDAR_INVVGLUUM39'}
%     {'NDAR_INVVD3KF4HE'}
%     {'NDAR_INVLW127JFT'}
%     {'NDAR_INVLRHY4MJ7'}
%     {'NDAR_INVU5DCN232'}
%     {'NDAR_INVYWEY5YL9'}
%     {'NDAR_INV56KML05T'}
%     {'NDAR_INVKVMU7J2K'}
%     {'NDAR_INVD32LCBUY'}
%     {'NDAR_INVJP0XXWG5'}
%     {'NDAR_INVRLH4F4JA'}
%     {'NDAR_INVZEV3PCPL'}
%     {'NDAR_INVR9Z48Y9H'}
%     {'NDAR_INVZ2E8YJN2'}
%     {'NDAR_INVZ3XXUKWN'}
%     {'NDAR_INVJ67B8F0J'}
%     {'NDAR_INVLBELG3A6'}
%     {'NDAR_INVNDEK4K38'}
%% 29 subjects were reset their sibling type and family ID 
%% according to the kinship in step 4 starting from code line 357
% LPC2, 3BZ3, M3Y9 --> 'NS'
% 88CW and X4DY, 7LW2 and 0W26, 5895 and 4CGV, ML75 and COU4, 85TU and BEZ6, R9C4 and 9HZA, KCYL and G6BN --> 'FS'
% HNGB and WY31 --> 'DZ' (ages were the same)
% RV36 and YJ7F, 7T9X and KA9M, HRKP and 6461, RD6E and LPZT, D51D and X4EW --> 'HS'
% the famID of the first subject in each pair was reset to use the famID of
% the second subject in this pair.
%% tabulate(imagingInfo2.sibtype)
%   Value    Count   Percent
%       0     2695     77.67%
%      10       17      0.49%
%     100      167      4.81%
%    1000      351     10.12%
%   10000      240      6.92%
% tabulate(imagingInfo2.famtype)
%   Value    Count   Percent
%       0     2694     87.55%
%      20        7      0.23%
%      30        1      0.03%
%     200       82      2.66%
%    2000      167      5.43%
%    2100        3      0.10%
%    3000        2      0.06%
%    4000        1      0.03%
%   20000      118      3.83%
%   20100        1      0.03%
%   21000        1      0.03%


%%%%%%%%%%%%
%% Step 1: load the questionnaire info for family relatedness
%%%%%%%%%%%%

fam = readtable('acspsw02.txt'); %family structure data,release 1.1
var = table2array(fam(1,:))';
ID = table2array(fam(2:end,4));

%3914 unique families
famID = str2double(table2array(fam(2:end,31)));
disp(size(unique(famID)));
%3types:1, 2, 3 
groupID = str2double(table2array(fam(2:end,32)));
%0=single,1=sibling,2=twin,3=triplet
famRela = str2double(table2array(fam(2:end,33)));
sexSame = str2double(table2array(fam(2:end,34)));
age = str2double(table2array(fam(2:end,7)));


%%%%%%%%%%%%
%% Step 2: load the kinship for the subjects selected in our study
%%%%%%%%%%%%
%% kinship values for each subject
kinship0884  = readtable('kinship_0.0884.txt');

imagingInfo = load('familyinfo.mat');
tic;
withinIndex = zeros(size(kinship0884,1),1);
for j = 1 : size(kinship0884,1)
    if mod(j,100) == 0
        disp(j);
    end
    if sum(strcmp(imagingInfo.c,kinship0884.ID1{j})) && sum(strcmp(imagingInfo.c,kinship0884.ID2{j}))
       withinIndex(j) = 1;
    end 
end
disp(toc);

kinship0884_within = kinship0884(withinIndex==1,:);
tic;
clear kinship
for i = 1 : length(imagingInfo.c)
    if mod(i,100) == 0
        disp(i);
    end
    kinship{i} = [];  
    idkin = find(strcmp(kinship0884_within.ID1, imagingInfo.c{i}) + strcmp(kinship0884_within.ID2, imagingInfo.c{i})); 
    kinship{i} = [kinship{i}, kinship0884_within.Kinship(idkin)];
end
disp(toc);

save('kinshipWithin.mat','kinship')

numsib = zeros(length(imagingInfo.c),1);
maxkin = zeros(length(imagingInfo.c),1);
for i = 1 : length(imagingInfo.c)
    numsib(i) = length(kinship{i});
    if numsib(i) > 0
        maxkin(i) = max(kinship{i});
    end
end
tabulate(numsib)
% >> tabulate(numsib)
%   Value    Count   Percent
%       0     2732     77.72%
%       1      758     21.56%
%       2       21      0.60%
%       3        4      0.11%

% threerecords = find(numsib==3);
% disp(imagingInfo.c(threerecords));
% % 
% %     {'NDAR_INVEC2CN745'}
% %     {'NDAR_INVGUAT6C5H'}
% %     {'NDAR_INVV9BM1RB5'}
% %     {'NDAR_INVVR0KF4XT'}
% % this family has two pairs of DZ twins in the imageing dataset 
% disp(famID(find(strcmp(ID, 'NDAR_INVEC2CN745'))));
% % 1415 is the family ID 
% clear idnow
% for i = 1 : length(threerecords)
%     idnow(i) = find(strcmp(ID, imagingInfo.c(threerecords)));
% end
% table(ID(idnow), age(idnow), famID(idnow), 'VariableNames', {'ID', 'age', 'fID'})

%%%%%%%%%%%%
%% Step 3: count the mismatches
%%%%%%%%%%%%

% combine the questionnair information and the genetic kinship 
%famRela: 0=single,1=sibling,2=twin,3=triplet
genetics  = readtable('ABCD_release_2.0.1.csv');
genemissing = zeros(length(imagingInfo.c),1);
count = zeros(1,4);

clear sibtypeName
record_all = [];
for i = 1 : length(imagingInfo.c)
    idnow  =  find(strcmp(ID, imagingInfo.c{i}));
    if maxkin(i) >= 0.354
        % mz twin
        sibtypeName(i,:) = 'MZ';
        if famRela(idnow) <= 1
            count(1) = count(1) + 1;
            disp(['inconsistent: MZ twin ', imagingInfo.c{i}, '; famID: ', num2str(famID(idnow))]);           
            idkin = find(strcmp(kinship0884_within.ID1, imagingInfo.c{i}) + strcmp(kinship0884_within.ID2, imagingInfo.c{i}));
            if ~isempty(setdiff(idkin, record_all))
                record_all = [record_all; idkin];
            end 
        end
    elseif maxkin(i) >= 0.177
        % dz twin or full sibling
        if famRela(idnow) > 1
            % dz twin            
            sibtypeName(i,:) = 'DZ';
        elseif famRela(idnow) == 1
            % full sibling
            sibtypeName(i,:) = 'FS';
        else
            count(2) = count(2) + 1;
            disp(['inconsistent: DZ or FS ', imagingInfo.c{i}, '; famID: ',  num2str(famID(idnow))]);
            
            idkin = find(strcmp(kinship0884_within.ID1, imagingInfo.c{i}) + strcmp(kinship0884_within.ID2, imagingInfo.c{i}));
            if ~isempty(setdiff(idkin, record_all))
                record_all = [record_all; idkin];
            end
        end
    elseif maxkin(i) >= 0.0884
        % half  sibling
        if famRela(idnow) == 1
            sibtypeName(i,:) = 'HS';
        else
            count(3) = count(3) + 1;
            disp(['inconsistent: HS ', imagingInfo.c{i}, '; famID: ',  num2str(famID(idnow))]);
            idkin = find(strcmp(kinship0884_within.ID1, imagingInfo.c{i}) + strcmp(kinship0884_within.ID2, imagingInfo.c{i}));            
            if ~isempty(setdiff(idkin, record_all))
                record_all = [record_all; idkin];
            end
        end
    else
        if famRela(idnow) == 0
            sibtypeName(i,:) = 'NS';
        else
            if sum(imagingInfo.famID_1==imagingInfo.famID_1(i)) == 1
                % if this subject was the only one selected from the family
                sibtypeName(i,:) = 'NS';
            else                      
                if famRela(idnow) == 1
                    sibtypeName(i,:) = 'SB';                    
                elseif famRela(idnow) == 2
                    sibtypeName(i,:) = 'TW';
                elseif famRela(idnow) == 3
                    sibtypeName(i,:) = '3W';
                end
                count(4) = count(4) + 1;
                disp(['inconsistent: NS ', imagingInfo.c{i}, '; famID: ',  num2str(famID(idnow))]);
                % the gene data are missing
                if ~ismember(imagingInfo.c{i},genetics.Var2)
                    genemissing(i) = 1; % the subject did not have gene data 
                else
                    genemissing(i) = 2; % the subject had gene data but the relatives might or might not                   
                end
            end
        end       
    end
end


record_all  = unique(record_all);
record_fam = [];
for ii = 1 : length(record_all)
    clear id1 id2
    id1 = find(strcmp(ID,  kinship0884_within.ID1(record_all(ii))));
    record_fam = [record_fam, id1];
    id2 = find(strcmp(ID,  kinship0884_within.ID2(record_all(ii))));
    record_fam = [record_fam, id2];
end
record_fam = unique(record_fam);
disp('--------------------------');
disp(kinship0884_within(record_all,:));
disp('--------------------------');
disp(fam(record_fam+1,[4:8,31:35]));
disp('--------------------------');    


%%
%        FID1               ID1               FID2               ID2              N_SNP       HetHet     IBS0     Kinship
%     ___________    __________________    ___________    __________________    __________    ______    ______    _______
% 
%     'AB0000239'    'NDAR_INVUDUKRV36'    'AB0000150'    'NDAR_INV87WRYJ7F'    5.0335e+05    0.0977    0.0262     0.096 
%     'AB0004657'    'NDAR_INVGBLP7T9X'    'AB0001790'    'NDAR_INV0RHLKA9M'    5.0262e+05    0.1084    0.0226    0.1229 
%     'AB0003273'    'NDAR_INV3ZM8HRKP'    'AB0003288'    'NDAR_INV93AY6461'    4.9907e+05    0.1007    0.0206    0.1266 
%     'AB0000291'    'NDAR_INVN0BDRD6E'    'AB0000292'    'NDAR_INVPM21LPZT'    5.0447e+05    0.1073    0.0219    0.1274 
%     'AB0003974'    'NDAR_INV23XZD51D'    'AB0003973'    'NDAR_INVBDMDX4EW'    4.9831e+05    0.1078    0.0191    0.1345 
%     'AB0004041'    'NDAR_INV6R35KCYL'    'AB0004394'    'NDAR_INV22E7G6BN'    5.0022e+05    0.1397    0.0148    0.2315 
%     'AB0005006'    'NDAR_INV08K0R9C4'    'AB0005008'    'NDAR_INVEJ0L9HZA'     4.874e+05    0.1396    0.0107     0.245 
%     'AB0002686'    'NDAR_INV6YFHHNGB'    'AB0002685'    'NDAR_INVGALZWY31'    4.9332e+05    0.1412    0.0096    0.2477 
%     'AB0001671'    'NDAR_INV1K1285TU'    'AB0002849'    'NDAR_INV1E35BEZ6'    5.0353e+05    0.1423    0.0098    0.2632 
%     'AB0001519'    'NDAR_INVJVWCML75'    'AB0003899'    'NDAR_INVK0LXC0U4'    4.9857e+05    0.1456    0.0101     0.269 
%     'AB0004407'    'NDAR_INV3YRA5895'    'AB0000828'    'NDAR_INV7LB64CGV'    4.9662e+05    0.1481    0.0094    0.2711 
%     'AB0004858'    'NDAR_INVFRT17LW2'    'AB0004861'    'NDAR_INVWZCU0W26'    5.0487e+05    0.1462    0.0091    0.2765 
%     'AB0002552'    'NDAR_INVK1CF88CW'    'AB0004370'    'NDAR_INVJZ6BX4DY'    5.0124e+05    0.1468    0.0085    0.2779 
%     'AB0003219'    'NDAR_INVTHEX3U3X'    'AB0003220'    'NDAR_INVBGMH3BZ3'    4.9616e+05    0.2284         0    0.4898 
%     'AB0005019'    'NDAR_INVF0GZLPC2'    'AB0018489'    'NDAR_INVGT92D3H6'    5.0452e+05    0.2306         0    0.4988

% count =
%     3    16    10    81 

%%%%%%%%%%%%
%% Step 4: deal with each category of mismatches in details
%%%%%%%%%%%%


%% deal with count(1) --> 3 : kinship was higher than 0.354 but family info was not twin or triplets
%         subjectkey          src_subject_id      interview_date    interview_age    gender    rel_family_id    rel_group_id    rel_relationship    rel_same_sex    race_ethnicity
%     __________________    __________________    ______________    _____________    ______    _____________    ____________    ________________    ____________    ______________
%
%     'NDAR_INVF0GZLPC2'    'NDAR_INVF0GZLPC2'     '08/22/2017'         '128'         'F'         '1583'            '1'               '0'               ''               '1'      
%     'NDAR_INVGT92D3H6'    'NDAR_INVGT92D3H6'     '08/22/2017'         '124'         'F'         '1595'            '1'               '0'               ''               '1'      
%     'NDAR_INVTHEX3U3X'    'NDAR_INVTHEX3U3X'     '06/13/2017'         '123'         'M'         '1374'            '1'               '2'               '1'              '1'      
%     'NDAR_INVNNA7JH41'    'NDAR_INVNNA7JH41'     '06/13/2017'         '123'         'M'         '1374'            '1'               '2'               '1'              '1'   
%     'NDAR_INVBGMH3BZ3'    'NDAR_INVBGMH3BZ3'     '06/09/2017'         '108'         'F'         '1265'            '1'               '0'               ''               '1'      

% LPC2 and D3H6 shared a kinship as 0.4988, but they were at different ages (128 and 124, respectively) and with different familyID
%% delete D3H6 and set sibtype(1672) = 'NS';    
deleteSet = find(strcmp(imagingInfo.c, 'NDAR_INVGT92D3H6'));  
sibtypeName(find(strcmp(imagingInfo.c, 'NDAR_INVF0GZLPC2')),:) = 'NS';
% 3U3X and 3BZ3 shared a kinship as 0.4898, but they were at different ages (123 and 108, respectively) and with different familyID


%% delete 3U3X and JH41 (according to the questionnaire they were twins, but JH41 did not have genetic data, so we could not determin MZ and DZ)
deleteSet = [deleteSet, find(strcmp(imagingInfo.c, 'NDAR_INVTHEX3U3X')), find(strcmp(imagingInfo.c, 'NDAR_INVNNA7JH41'))]; %find(strcmp(genetics.Var2, 'NDAR_INVNNA7JH41'))
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INVBGMH3BZ3')==1,:) = 'NS';

%% deal with count(2) --> 16 : kinship was higher than 0.177 but family info was not twin or triplets or sibling
% if ages were the same, then DZ; otherwise, FS;
% the family ID was reset to the same family ID for each pair
%         subjectkey          src_subject_id      interview_date    interview_age    gender    rel_family_id    rel_group_id    rel_relationship    rel_same_sex    race_ethnicity
%     __________________    __________________    ______________    _____________    ______    _____________    ____________    ________________    ____________    ______________
% 
%     'NDAR_INVK1CF88CW'    'NDAR_INVK1CF88CW'     '05/03/2017'         '131'         'M'         '3130'            '1'               '0'               ''               '1'      
%     'NDAR_INVJZ6BX4DY'    'NDAR_INVJZ6BX4DY'     '07/29/2017'         '112'         'F'         '3129'            '1'               '0'               ''               '1'      
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INVK1CF88CW')==1,:) = 'FS';
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INVJZ6BX4DY')==1,:) = 'FS';
imagingInfo.famID_1(strcmp(imagingInfo.c, 'NDAR_INVK1CF88CW')==1) = imagingInfo.famID_1(strcmp(imagingInfo.c, 'NDAR_INVJZ6BX4DY')==1); 

%     'NDAR_INVFRT17LW2'    'NDAR_INVFRT17LW2'     '08/20/2017'         '110'         'M'         '3805'            '1'               '0'               ''               '1'      
%     'NDAR_INVWZCU0W26'    'NDAR_INVWZCU0W26'     '08/20/2017'         '129'         'F'         '3855'            '1'               '0'               ''               '1'      
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INVFRT17LW2')==1,:) = 'FS';
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INVWZCU0W26')==1,:) = 'FS';
imagingInfo.famID_1(strcmp(imagingInfo.c, 'NDAR_INVFRT17LW2')==1) = imagingInfo.famID_1(strcmp(imagingInfo.c, 'NDAR_INVWZCU0W26')==1); 

%     'NDAR_INV3YRA5895'    'NDAR_INV3YRA5895'     '07/31/2017'         '109'         'F'         '2988'            '1'               '0'               ''               '1'      
%     'NDAR_INV7LB64CGV'    'NDAR_INV7LB64CGV'     '01/05/2017'         '130'         'F'         '3012'            '1'               '0'               ''               '1'      
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INV3YRA5895')==1,:) = 'FS';
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INV7LB64CGV')==1,:) = 'FS';
imagingInfo.famID_1(strcmp(imagingInfo.c, 'NDAR_INV3YRA5895')==1) = imagingInfo.famID_1(strcmp(imagingInfo.c, 'NDAR_INV7LB64CGV')==1); 

%     'NDAR_INVJVWCML75'    'NDAR_INVJVWCML75'     '02/26/2017'         '130'         'M'         '4443'            '1'               '0'               ''               '1'      
%     'NDAR_INVK0LXC0U4'    'NDAR_INVK0LXC0U4'     '07/10/2017'         '108'         'M'         '4445'            '1'               '0'               ''               '1'      
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INVJVWCML75')==1,:) = 'FS';
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INVK0LXC0U4')==1,:) = 'FS';
imagingInfo.famID_1(strcmp(imagingInfo.c, 'NDAR_INVJVWCML75')==1) = imagingInfo.famID_1(strcmp(imagingInfo.c, 'NDAR_INVK0LXC0U4')==1); 

%     'NDAR_INV1K1285TU'    'NDAR_INV1K1285TU'     '02/25/2017'         '128'         'M'         '1459'            '1'               '0'               ''               '1'      
%     'NDAR_INV1E35BEZ6'    'NDAR_INV1E35BEZ6'     '05/20/2017'         '109'         'F'         '1456'            '1'               '0'               ''               '1'      
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INV1K1285TU')==1,:) = 'FS';
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INV1E35BEZ6')==1,:) = 'FS';
imagingInfo.famID_1(strcmp(imagingInfo.c, 'NDAR_INV1K1285TU')==1) = imagingInfo.famID_1(strcmp(imagingInfo.c, 'NDAR_INV1E35BEZ6')==1); 

%     'NDAR_INV6YFHHNGB'    'NDAR_INV6YFHHNGB'     '05/13/2017'         '125'         'F'         '4177'            '1'               '0'               ''               '2'      
%     'NDAR_INVGALZWY31'    'NDAR_INVGALZWY31'     '05/13/2017'         '125'         'F'         '4225'            '1'               '0'               ''               '2'      
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INV6YFHHNGB')==1,:) = 'DZ';
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INVGALZWY31')==1,:) = 'DZ';
imagingInfo.famID_1(strcmp(imagingInfo.c, 'NDAR_INV6YFHHNGB')==1) = imagingInfo.famID_1(strcmp(imagingInfo.c, 'NDAR_INVGALZWY31')==1); 

%     'NDAR_INV08K0R9C4'    'NDAR_INV08K0R9C4'     '08/23/2017'         '127'         'M'         '3266'            '1'               '0'               ''               '1'   
%     'NDAR_INVEJ0L9HZA'    'NDAR_INVEJ0L9HZA'     '08/23/2017'         '109'         'M'         '3320'            '1'               '0'               ''               '1'      
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INV08K0R9C4')==1,:) = 'FS';
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INVEJ0L9HZA')==1,:) = 'FS';
imagingInfo.famID_1(strcmp(imagingInfo.c, 'NDAR_INV08K0R9C4')==1) = imagingInfo.famID_1(strcmp(imagingInfo.c, 'NDAR_INVEJ0L9HZA')==1); 

%     'NDAR_INV6R35KCYL'    'NDAR_INV6R35KCYL'     '07/14/2017'         '127'         'F'         '3615'            '1'               '0'               ''               '3'      
%     'NDAR_INV22E7G6BN'    'NDAR_INV22E7G6BN'     '07/14/2017'         '116'         'F'         '3622'            '1'               '0'               ''               '3'      
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INV6R35KCYL')==1,:) = 'FS';
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INV22E7G6BN')==1,:) = 'FS';
imagingInfo.famID_1(strcmp(imagingInfo.c, 'NDAR_INV6R35KCYL')==1) = imagingInfo.famID_1(strcmp(imagingInfo.c, 'NDAR_INV22E7G6BN')==1); 


%% deal with count(3) --> 10 : kinship was higher than 0.0884 but family info was not sibling

%        FID1               ID1               FID2               ID2              N_SNP       HetHet     IBS0     Kinship
%     ___________    __________________    ___________    __________________    __________    ______    ______    _______
% 
%     'AB0000239'    'NDAR_INVUDUKRV36'    'AB0000150'    'NDAR_INV87WRYJ7F'    5.0335e+05    0.0977    0.0262     0.096 
%     'AB0004657'    'NDAR_INVGBLP7T9X'    'AB0001790'    'NDAR_INV0RHLKA9M'    5.0262e+05    0.1084    0.0226    0.1229 
%     'AB0003273'    'NDAR_INV3ZM8HRKP'    'AB0003288'    'NDAR_INV93AY6461'    4.9907e+05    0.1007    0.0206    0.1266 
%     'AB0000291'    'NDAR_INVN0BDRD6E'    'AB0000292'    'NDAR_INVPM21LPZT'    5.0447e+05    0.1073    0.0219    0.1274 
%     'AB0003974'    'NDAR_INV23XZD51D'    'AB0003973'    'NDAR_INVBDMDX4EW'    4.9831e+05    0.1078    0.0191    0.1345 

%     'NDAR_INVUDUKRV36'    'NDAR_INVUDUKRV36'     '10/25/2016'         '127'         'M'         '630'             '1'               '0'               ''               '1'      
%     'NDAR_INV87WRYJ7F'    'NDAR_INV87WRYJ7F'     '10/06/2016'         '129'         'M'         '392'             '1'               '0'               ''               '1'      
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INVUDUKRV36')==1,:) = 'HS';
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INV87WRYJ7F')==1,:) = 'HS';
imagingInfo.famID_1(strcmp(imagingInfo.c, 'NDAR_INVUDUKRV36')==1) = imagingInfo.famID_1(strcmp(imagingInfo.c, 'NDAR_INV87WRYJ7F')==1); 

%     'NDAR_INVGBLP7T9X'    'NDAR_INVGBLP7T9X'     '08/10/2017'         '110'         'F'         '2437'            '1'               '0'               ''               '3'      
%     'NDAR_INV0RHLKA9M'    'NDAR_INV0RHLKA9M'     '03/03/2017'         '129'         'F'         '2278'            '1'               '0'               ''               '3'      
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INVGBLP7T9X')==1,:) = 'HS';
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INV0RHLKA9M')==1,:) = 'HS';
imagingInfo.famID_1(strcmp(imagingInfo.c, 'NDAR_INVGBLP7T9X')==1) = imagingInfo.famID_1(strcmp(imagingInfo.c, 'NDAR_INV0RHLKA9M')==1); 

%     'NDAR_INV3ZM8HRKP'    'NDAR_INV3ZM8HRKP'     '06/14/2017'         '124'         'M'         '3738'            '1'               '0'               ''               '1'      
%     'NDAR_INV93AY6461'    'NDAR_INV93AY6461'     '06/15/2017'         '131'         'F'         '3770'            '1'               '0'               ''               '1'      
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INV3ZM8HRKP')==1,:) = 'HS';
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INV93AY6461')==1,:) = 'HS';
imagingInfo.famID_1(strcmp(imagingInfo.c, 'NDAR_INV3ZM8HRKP')==1) = imagingInfo.famID_1(strcmp(imagingInfo.c, 'NDAR_INV93AY6461')==1); 

%     'NDAR_INVN0BDRD6E'    'NDAR_INVN0BDRD6E'     '11/05/2016'         '109'         'F'         '3162'            '1'               '0'               ''               '2'      
%     'NDAR_INVPM21LPZT'    'NDAR_INVPM21LPZT'     '11/05/2016'         '131'         'M'         '3175'            '1'               '0'               ''               '2'      
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INVN0BDRD6E')==1,:) = 'HS';
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INVPM21LPZT')==1,:) = 'HS';
imagingInfo.famID_1(strcmp(imagingInfo.c, 'NDAR_INVN0BDRD6E')==1) = imagingInfo.famID_1(strcmp(imagingInfo.c, 'NDAR_INVPM21LPZT')==1); 

%     'NDAR_INV23XZD51D'    'NDAR_INV23XZD51D'     '07/15/2017'         '131'         'M'         '2089'            '1'               '0'               ''               '2'      
%     'NDAR_INVBDMDX4EW'    'NDAR_INVBDMDX4EW'     '07/15/2017'         '118'         'M'         '2131'            '1'               '0'               ''               '2'      
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INV23XZD51D')==1,:) = 'HS';
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INVBDMDX4EW')==1,:) = 'HS';
imagingInfo.famID_1(strcmp(imagingInfo.c, 'NDAR_INV23XZD51D')==1) = imagingInfo.famID_1(strcmp(imagingInfo.c, 'NDAR_INVBDMDX4EW')==1); 

%% 81 subjects had siblings or twins as specified in the questionnaire but did not have the kinship estimation owing to either
%% no gene data of their owns (40) or no gene data of their relatives (41)
%% among them 34--registerred with siblins, 47--registerred with twins
% >> tabulate(sibtypeName)
%   Value    Count   Percent
%      FS      154      4.41%
%      NS     2651     75.98%
%      DZ      351     10.06%
%      MZ      244      6.99%
%      TW       47      1.35%
%      HS        8      0.23%
%      SB       34      0.97%
sibtypeNameCell = mat2cell(sibtypeName,ones(length(sibtypeName),1),2);

%% only one from each family remain: the one with rel_group_id == 2 was deleted
sbid = find(strcmp(sibtypeNameCell, 'SB'));
idx = [];
for i = 1 : length(sbid)
    idx = [idx, find(strcmp(ID, imagingInfo.c(sbid(i))))];
    if (groupID(idx(i)) == 2)        
        deleteSet = [deleteSet, sbid(i)];
    end
end
disp(sortrows(fam(idx+1,[4:8,31:35]), 'rel_family_id'));
sibtypeName(sbid,:) = repmat('NS',length(sbid),1);

%%
%         subjectkey          src_subject_id      interview_date    interview_age    gender    rel_family_id    rel_group_id    rel_relationship    rel_same_sex    race_ethnicity
%     __________________    __________________    ______________    _____________    ______    _____________    ____________    ________________    ____________    ______________
% 
%     'NDAR_INV96WYY4C5'    'NDAR_INV96WYY4C5'     '06/21/2017'         '121'         'F'         '2107'            '1'               '1'                ''              '1'      
%     'NDAR_INVYDEWFLED'    'NDAR_INVYDEWFLED'     '06/21/2017'         '108'         'M'         '2107'            '2'               '1'                ''              '1'      
%     'NDAR_INV6NRB61V1'    'NDAR_INV6NRB61V1'     '05/30/2017'         '109'         'M'         '2210'            '2'               '1'                ''              '2'      
%     'NDAR_INVMZKNKER7'    'NDAR_INVMZKNKER7'     '05/30/2017'         '125'         'F'         '2210'            '1'               '1'                ''              '2'      
%     'NDAR_INV75H3LM5P'    'NDAR_INV75H3LM5P'     '07/09/2017'         '111'         'M'         '3046'            '1'               '1'                ''              '2'      
%     'NDAR_INVDLMML76V'    'NDAR_INVDLMML76V'     '07/09/2017'         '111'         'M'         '3046'            '2'               '1'                ''              '2'      
%     'NDAR_INVB3438R23'    'NDAR_INVB3438R23'     '04/13/2017'         '122'         'F'         '3070'            '2'               '1'                ''              '1'      
%     'NDAR_INVB8234KWV'    'NDAR_INVB8234KWV'     '07/03/2017'         '111'         'M'         '3070'            '1'               '1'                ''              '4'      
%     'NDAR_INVF19VWPA1'    'NDAR_INVF19VWPA1'     '12/01/2016'         '129'         'F'         '3123'            '2'               '1'                ''              '5'      
%     'NDAR_INVHD40HGZ7'    'NDAR_INVHD40HGZ7'     '12/01/2016'         '108'         'M'         '3123'            '1'               '1'                ''              '5'      
%     'NDAR_INVVB04BHVZ'    'NDAR_INVVB04BHVZ'     '01/16/2017'         '131'         'M'         '3212'            '1'               '1'                ''              '5'      
%     'NDAR_INVZFYGCHKB'    'NDAR_INVZFYGCHKB'     '01/16/2017'         '108'         'M'         '3212'            '2'               '1'                ''              '5'      
%     'NDAR_INVE59K1A64'    'NDAR_INVE59K1A64'     '08/26/2017'         '109'         'M'         '3317'            '1'               '1'                ''              '1'      
%     'NDAR_INVVF539UBP'    'NDAR_INVVF539UBP'     '07/20/2017'         '131'         'F'         '3317'            '2'               '1'                ''              '3'      
%     'NDAR_INV31JMTLE9'    'NDAR_INV31JMTLE9'     '08/17/2017'         '115'         'M'         '3371'            '2'               '1'                ''              '2'      
%     'NDAR_INVX46YGYN7'    'NDAR_INVX46YGYN7'     '08/17/2017'         '121'         'M'         '3371'            '1'               '1'                ''              '2'      
%     'NDAR_INV4GHXF68C'    'NDAR_INV4GHXF68C'     '04/17/2017'         '109'         'F'         '3403'            '1'               '1'                ''              '4'      
%     'NDAR_INVAEMLZKT4'    'NDAR_INVAEMLZKT4'     '04/17/2017'         '114'         'F'         '3403'            '2'               '1'                ''              '4'      
%     'NDAR_INV7CN47GEF'    'NDAR_INV7CN47GEF'     '01/21/2017'         '123'         'M'         '4004'            '2'               '1'                ''              '3'      
%     'NDAR_INVCGWV06JE'    'NDAR_INVCGWV06JE'     '01/21/2017'         '109'         'F'         '4004'            '1'               '1'                ''              '3'      
%     'NDAR_INVJCE41W2C'    'NDAR_INVJCE41W2C'     '07/01/2017'         '110'         'M'         '4233'            '1'               '1'                ''              '1'      
%     'NDAR_INVLDH3YU2R'    'NDAR_INVLDH3YU2R'     '07/01/2017'         '119'         'F'         '4233'            '2'               '1'                ''              '1'      
%     'NDAR_INVMWC4BV86'    'NDAR_INVMWC4BV86'     '11/17/2016'         '131'         'M'         '571'             '2'               '1'                ''              '1'      
%     'NDAR_INVN5FRRMAK'    'NDAR_INVN5FRRMAK'     '07/21/2017'         '110'         'F'         '571'             '1'               '1'                ''              '1'      
%     'NDAR_INVNT6BJVLB'    'NDAR_INVNT6BJVLB'     '04/28/2017'         '127'         'M'         '574'             '1'               '1'                ''              '2'      
%     'NDAR_INVWDCENDZ9'    'NDAR_INVWDCENDZ9'     '04/28/2017'         '123'         'F'         '574'             '2'               '1'                ''              '1'      
%     'NDAR_INVAUR9UZH4'    'NDAR_INVAUR9UZH4'     '11/03/2016'         '128'         'F'         '643'             '2'               '1'                ''              '5'      
%     'NDAR_INVV9D0WWZM'    'NDAR_INVV9D0WWZM'     '11/03/2016'         '114'         'M'         '643'             '1'               '1'                ''              '5'      
%     'NDAR_INVC6MPJ7ZF'    'NDAR_INVC6MPJ7ZF'     '01/26/2017'         '131'         'F'         '663'             '2'               '1'                ''              '1'      
%     'NDAR_INVWXLHTR31'    'NDAR_INVWXLHTR31'     '08/17/2017'         '120'         'F'         '663'             '1'               '1'                ''              '1'      
%     'NDAR_INV62EUKL2P'    'NDAR_INV62EUKL2P'     '05/17/2017'         '108'         'M'         '917'             '1'               '1'                ''              '1'      
%     'NDAR_INVWMM5T2PL'    'NDAR_INVWMM5T2PL'     '05/17/2017'         '131'         'F'         '917'             '2'               '1'                ''              '1'      
%     'NDAR_INVE7CFEXUC'    'NDAR_INVE7CFEXUC'     '01/25/2017'         '110'         'F'         '92'              '2'               '1'                ''              '1'      
%     'NDAR_INVGFLA5X8W'    'NDAR_INVGFLA5X8W'     '02/08/2017'         '125'         'M'         '92'              '1'               '1'                ''              '1'

%% only one from each family remain: the one with rel_group_id == 2 was deleted
twid = find(strcmp(sibtypeNameCell, 'TW'));
idx = [];
for i = 1 : length(twid)
    idx = [idx, find(strcmp(ID, imagingInfo.c(twid(i))))];
end
disp(sortrows(fam(idx+1,[4:8,31:35]), 'rel_family_id'));

fidnow = unique(famID(idx));
for i = 1 : length(fidnow)
    setnow = find(famID(idx)==fidnow(i));
    if length(setnow) == 2
        deleteSet = [deleteSet, twid(setnow(2))]; % delete one of the twin whose kinship could not be estimated owing to missing genetic data 
    end
end
sibtypeName(twid,:) = repmat('NS', length(twid), 1);




%        subjectkey          src_subject_id      interview_date    interview_age    gender    rel_family_id    rel_group_id    rel_relationship    rel_same_sex    race_ethnicity
%     __________________    __________________    ______________    _____________    ______    _____________    ____________    ________________    ____________    ______________
% 
%     'NDAR_INV5L4LYVBY'    'NDAR_INV5L4LYVBY'     '03/26/2017'         '129'         'M'         '1194'            '1'               '2'               '1'              '5'      
%     'NDAR_INVEKN9FWNN'    'NDAR_INVEKN9FWNN'     '03/26/2017'         '129'         'M'         '1194'            '1'               '2'               '1'              '5'      
%     'NDAR_INV8MMU6YZ0'    'NDAR_INV8MMU6YZ0'     '04/06/2017'         '130'         'F'         '1253'            '1'               '2'               '1'              '1'      
%     'NDAR_INVA6FR087J'    'NDAR_INVA6FR087J'     '04/06/2017'         '130'         'F'         '1253'            '1'               '2'               '1'              '1'      
%     'NDAR_INV8C5PMGAM'    'NDAR_INV8C5PMGAM'     '01/21/2017'         '119'         'M'         '1277'            '1'               '2'               '1'              '2'      
%     'NDAR_INVCTH2AK0B'    'NDAR_INVCTH2AK0B'     '01/21/2017'         '119'         'M'         '1277'            '1'               '2'               '1'              '2'      
%     'NDAR_INVFD93F59C'    'NDAR_INVFD93F59C'     '07/06/2017'         '118'         'M'         '1317'            '1'               '2'               '1'              '2'      
%     'NDAR_INVTWU5RFVX'    'NDAR_INVTWU5RFVX'     '07/06/2017'         '118'         'M'         '1317'            '1'               '2'               '1'              '2'      
%     'NDAR_INV8HX2TBDP'    'NDAR_INV8HX2TBDP'     '08/21/2017'         '122'         'M'         '1330'            '1'               '2'               '1'              '1'      
%     'NDAR_INVH0DGUW54'    'NDAR_INVH0DGUW54'     '08/21/2017'         '122'         'M'         '1330'            '1'               '2'               '1'              '1'      
%     'NDAR_INVPDV2UAGX'    'NDAR_INVPDV2UAGX'     '07/26/2017'         '115'         'F'         '1379'            '1'               '2'               '1'              '5'      
%     'NDAR_INVVGLUUM39'    'NDAR_INVVGLUUM39'     '07/26/2017'         '115'         'F'         '1379'            '1'               '2'               '1'              '5'      
%     'NDAR_INV4MT8HH02'    'NDAR_INV4MT8HH02'     '07/12/2017'         '129'         'M'         '1416'            '1'               '2'               '1'              '1'      
%     'NDAR_INVVD3KF4HE'    'NDAR_INVVD3KF4HE'     '07/12/2017'         '129'         'M'         '1416'            '1'               '2'               '1'              '1'      
%     'NDAR_INV0PLKFP06'    'NDAR_INV0PLKFP06'     '02/10/2017'         '114'         'F'         '2020'            '1'               '2'               '1'              '1'      
%     'NDAR_INVLW127JFT'    'NDAR_INVLW127JFT'     '02/10/2017'         '114'         'F'         '2020'            '1'               '2'               '1'              '1'      
%     'NDAR_INVH8KHAFL7'    'NDAR_INVH8KHAFL7'     '06/30/2017'         '125'         'M'         '2199'            '1'               '2'               '1'              '1'      
%     'NDAR_INVLRHY4MJ7'    'NDAR_INVLRHY4MJ7'     '06/30/2017'         '125'         'M'         '2199'            '1'               '2'               '1'              '1'      
%     'NDAR_INVPG5HVR50'    'NDAR_INVPG5HVR50'     '01/23/2017'         '121'         'M'         '2222'            '1'               '2'               '1'              '1'      
%     'NDAR_INVU5DCN232'    'NDAR_INVU5DCN232'     '01/23/2017'         '121'         'M'         '2222'            '1'               '2'               '1'              '1'      
%     'NDAR_INV73J1L7TG'    'NDAR_INV73J1L7TG'     '03/16/2017'         '118'         'M'         '2588'            '1'               '2'               '0'              '4'      
%     'NDAR_INVYWEY5YL9'    'NDAR_INVYWEY5YL9'     '03/16/2017'         '118'         'F'         '2588'            '1'               '2'               '0'              '4'      
%     'NDAR_INVAZY0D7JP'    'NDAR_INVAZY0D7JP'     '07/03/2017'         '120'         'F'         '268'             '1'               '2'               '1'              '1'      
%     'NDAR_INVV7D70A79'    'NDAR_INVV7D70A79'     '07/03/2017'         '120'         'F'         '268'             '1'               '2'               '1'              '1'      
%     'NDAR_INV01EN91PG'    'NDAR_INV01EN91PG'     '04/19/2017'         '113'         'F'         '2699'            '1'               '2'               '1'              '1'      
%     'NDAR_INV56KML05T'    'NDAR_INV56KML05T'     '04/19/2017'         '113'         'F'         '2699'            '1'               '2'               '1'              '1'      
%     'NDAR_INV08EFDKZ6'    'NDAR_INV08EFDKZ6'     '01/27/2017'         '131'         'M'         '2703'            '1'               '2'               '1'              '2'      
%     'NDAR_INVKVMU7J2K'    'NDAR_INVKVMU7J2K'     '01/27/2017'         '131'         'M'         '2703'            '1'               '2'               '1'              '2'      
%     'NDAR_INV3RFU7ZRP'    'NDAR_INV3RFU7ZRP'     '07/10/2017'         '121'         'F'         '2820'            '1'               '2'               '1'              '1'      
%     'NDAR_INVD32LCBUY'    'NDAR_INVD32LCBUY'     '07/10/2017'         '121'         'F'         '2820'            '1'               '2'               '1'              '1'      
%     'NDAR_INVFPX9E85L'    'NDAR_INVFPX9E85L'     '12/01/2016'         '129'         'F'         '2871'            '1'               '2'               '1'              '5'      
%     'NDAR_INVJP0XXWG5'    'NDAR_INVJP0XXWG5'     '12/01/2016'         '129'         'F'         '2871'            '1'               '2'               '1'              '5'      
%     'NDAR_INVNX6KVGNJ'    'NDAR_INVNX6KVGNJ'     '08/12/2017'         '113'         'F'         '2912'            '1'               '2'               '1'              '1'      
%     'NDAR_INVRLH4F4JA'    'NDAR_INVRLH4F4JA'     '08/12/2017'         '113'         'F'         '2912'            '1'               '2'               '1'              '1'      
%     'NDAR_INV50411TBC'    'NDAR_INV50411TBC'     '06/23/2017'         '115'         'M'         '2951'            '1'               '2'               '1'              '1'      
%     'NDAR_INVZEV3PCPL'    'NDAR_INVZEV3PCPL'     '06/23/2017'         '115'         'M'         '2951'            '1'               '2'               '1'              '1'      
%     'NDAR_INV1285PMCK'    'NDAR_INV1285PMCK'     '05/31/2017'         '122'         'F'         '404'             '1'               '2'               '1'              '5'      
%     'NDAR_INVWCV58UKD'    'NDAR_INVWCV58UKD'     '05/09/2017'         '110'         'F'         '404'             '1'               '2'               '1'              '1'      
%     'NDAR_INV93U4YVJ1'    'NDAR_INV93U4YVJ1'     '07/06/2017'         '127'         'M'         '4480'            '1'               '2'               '1'              '1'      
%     'NDAR_INVR9Z48Y9H'    'NDAR_INVR9Z48Y9H'     '07/06/2017'         '127'         'M'         '4480'            '1'               '2'               '1'              '1'      
%     'NDAR_INVTXK16WY8'    'NDAR_INVTXK16WY8'     '07/02/2017'         '131'         'M'         '4485'            '1'               '2'               '1'              '1'      
%     'NDAR_INVZ2E8YJN2'    'NDAR_INVZ2E8YJN2'     '07/02/2017'         '131'         'M'         '4485'            '1'               '2'               '1'              '1'      
%     'NDAR_INV02EBX0JJ'    'NDAR_INV02EBX0JJ'     '12/01/2016'         '117'         'M'         '4515'            '1'               '2'               '1'              '1'      
%     'NDAR_INVZ3XXUKWN'    'NDAR_INVZ3XXUKWN'     '12/01/2016'         '117'         'M'         '4515'            '1'               '2'               '1'              '1'  

%%
%  'NDAR_INVJ67B8F0J'    'NDAR_INVJ67B8F0J'     '01/28/2017'         '129'         'F'         '1289'            '1'               '2'               '1'              '1'      
%  'NDAR_INVDTHJM3Y9'    'NDAR_INVDTHJM3Y9'     '01/28/2017'         '129'         'F'         '1289'            '1'               '2'               '1'              '1'      

%  'NDAR_INVNDEK4K38'    'NDAR_INVNDEK4K38'     '01/29/2017'         '123'         'M'         '1360'            '1'               '2'               '1'              '1'      
%  'NDAR_INVLBELG3A6'    'NDAR_INVLBELG3A6'     '01/29/2017'         '123'         'M'         '1360'            '1'               '2'               '1'              '1'      

%  AB0001110	NDAR_INVDTHJM3Y9	AB0001111	NDAR_INVNDEK4K38	493514	0.1401	0.0099	0.2474

deleteSet = [deleteSet, find(strcmp(imagingInfo.c,  'NDAR_INVJ67B8F0J')),  find(strcmp(imagingInfo.c,  'NDAR_INVLBELG3A6')), find(strcmp(imagingInfo.c,  'NDAR_INVNDEK4K38'))];
sibtypeName(strcmp(imagingInfo.c, 'NDAR_INVDTHJM3Y9')==1,:) = 'NS';

%               'NDAR_INVTHEX3U3X'    'NDAR_INVTHEX3U3X'     '06/13/2017'         '123'         'M'         '1374'            '1'               '2'               '1'              '1'      
%  [No  gene;  has been deleted]   'NDAR_INVNNA7JH41'    'NDAR_INVNNA7JH41'     '06/13/2017'         '123'         'M'         '1374'            '1'               '2'               '1'              '1'      
%%
imagingInfo2.c = imagingInfo.c(setdiff(1:length(imagingInfo.c),deleteSet),1);
imagingInfo2.famID_1 = imagingInfo.famID_1(setdiff(1:length(imagingInfo.c),deleteSet),1);
imagingInfo2.famRela = sibtypeName(setdiff(1:length(imagingInfo.c),deleteSet),:);

%% tabulate(imagingInfo2.famRela)
%   Value    Count   Percent
%      FS      167      4.82%
%      NS     2693     77.65%
%      DZ      351     10.12%
%      MZ      240      6.92%
%      HS       17      0.49%
%% family type
F = unique(imagingInfo2.famID_1);
% define a numeric sibling type
% assume that no family has more than 9 children
for i = 1 : length(imagingInfo2.c)
    switch imagingInfo2.famRela(i,:)
        case 'NS'
            imagingInfo2.sibtype(i) = 0;
        case 'HS'
            imagingInfo2.sibtype(i) = 10;
        case 'FS'
            imagingInfo2.sibtype(i) = 100;
        case 'DZ'
            imagingInfo2.sibtype(i) = 1000;
        case 'MZ'
            imagingInfo2.sibtype(i) = 10000;
    end          
end
% tabulate(imagingInfo2.sibtype)
%   Value    Count   Percent
%       0     2695     77.67%
%      10       17      0.49%
%     100      167      4.81%
%    1000      351     10.12%
%   10000      240      6.92%
%%  3077 families
for i = 1 : length(F)
    idx = find(imagingInfo2.famID_1 == F(i));
    imagingInfo2.famtype(i) = sum(imagingInfo2.sibtype(idx));
end
% tabulate(imagingInfo2.famtype)
%   Value    Count   Percent
%       0     2694     87.55%
%      20        7      0.23%
%      30        1      0.03%
%     200       82      2.66%
%    2000      167      5.43%
%    2100        3      0.10%
%    3000        2      0.06%
%    4000        1      0.03%
%   20000      118      3.83%
%   20100        1      0.03%
%   21000        1      0.03%
%% if twins are marked as 1st degree, MZ; otherwise DZ
save('imagingInfo2.mat', 'imagingInfo2');