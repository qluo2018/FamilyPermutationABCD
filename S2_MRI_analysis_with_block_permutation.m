%%%%%%%%%%%%%%
%% Step 2£ºBrain association analysis (using sleep as an example)
%%%%%%%%%%%%%%

%% Step 2.1 % Conduct permutation exchangable block.
% Code revised from hcp2blocks.m by Anderson M. Winkler
% http://brainder.org
% Reference:
% Winkler AM, Webster MA, Vidaurre D, Nichols TE, Smith SM.
% Multi-level block permutation. Neuroimage. 2015;123:253-68.

%% written by Shen Chun, cshen17@fudan.edu.cn
%% reviewed by Dr Qiang Luo, qluo@fudan.edu.cn
%% released on 21 Mar 2020
%% please cite: Shen, et al. Biological Psychiatry 2020

%load data
load('familyinfo.mat') % 3 variables: ID (N-by-1 cell), famRela(N-by-1 double), famID (N-by-1, double)

N = size(ID,1);

%sibtype
sibtype = zeros(N,1);
for n = 1:N,
    if famRela(n)<2,%single,sibling(not twin)=10
        sibtype(n) = 10;
    elseif famRela(n)==2,
        sibtype(n) = 100;%twin
    elseif famRela(n)==3,
        sibtype(n) = 1000;%triplet
    end
end

%sibtype based on genetic data
%sibtype = zeros(N,1);
%for n = 1:N,
%    if any(strcmpi(famRela(n,:),{'FS','NS'})),
%        sibtype(n) = 10;
%    elseif strcmpi(famRela(n,:),'HS'),
%        sibtype(n) = 11;
%    elseif strcmpi(famRela(n,:),'DZ'),
%        sibtype(n) = 100;
%    elseif strcmpi(famRela(n,:),'MZ'),
%        sibtype(n) = 1000;
%    end
%end

%famtype
Ufam = unique(famID);
famtype = zeros(N,1);
for f = 1:numel(Ufam),
    fidx = Ufam(f) == famID;
    famtype(fidx) = sum(sibtype(fidx));
end

%special subs, two pair of twins in one family
%no. 1607 and 3017 are twins
sibtype(1607)=101;
sibtype(3017)=101;

[~,idx] = sortrows([famID sibtype]);% ascending sort by famid(1) sibtype(2) age(3)
[~,idxback] = sort(idx);
sibtype = sibtype(idx);
famID   = famID(idx);
famtype = famtype(idx);

% Now make the blocks for each family
B = cell(numel(Ufam),1);
for f = 1:numel(Ufam),
    fidx = Ufam(f) == famID;
    ft = famtype(find(fidx,1));
    if any(ft == [10,20,200,400,3000]),
        B{f} = horzcat(famID(fidx),sibtype(fidx));
    %for genetic data
    %if any(ft == [10,20,22,200,300,400,2000]),% within-family shuffle
    %    B{f} = horzcat(famID(fidx),sibtype(fidx));
    else
        B{f} = horzcat(-famID(fidx),sibtype(fidx));
    end
end

% Concatenate all. Prepending the famtype ensures that the
% families of the same type can be shuffled whole-block. Also,
% add column with -1, for within-block at the outermost level
B = horzcat(-ones(N,1),famtype,cell2mat(B));

B = B(idxback,:);
            
save EB_abcd B ID;

%subject No. column
B(:,5) = 1:1:3515;

%create 5000 times permutation orders
[Pset, VG] = palm_quickperms([], B,5001);
Pset = Pset(:,2:end);
save permOrder Pset VG ID;
%% Step 2.2 Brain wide association with block permutation
% Whole-brain voxel-wise analysis
% Linear regression + multi-level block permutation-based cluster-level correction
% Permutation code was adapted from Weikang Gong (E-mail: weikanggong@gmail.com)

%load data
load('DATA.mat')
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

%parameters
CDT = 0.001;% cluster defining threshold p-value
X = GMV;
design = [dysomnia,cov_comb];
mask = mask_index;% mask
nperm = 5000;% nuber of permutation

[a1,a2,a3]=size(mask);
[nn,~]=size(X);
qq=size(design,2);
[d1,d2,d3]=ind2sub(size(mask),find(mask~=0));
dim=[d1,d2,d3];
mm=length(CDT);
thre1=abs(tinv(CDT,nn-qq-1));%positive t threshold
thre2=tinv(CDT,nn-qq-1);%negative t threshold
%Permutation SPM
maxcl1=zeros(nperm,mm);%maximum cluster size in each permutation
maxcl2=zeros(nperm,mm);

%permutation order considering family structure
load('permOrder.mat','Pset')
for i=1:nperm
    
    ind=Pset(:,i);%permutation order
    ts=BWAS_Tregression(design(ind,:),X);
    for j=1:mm
        ind1=ts>thre1(j);
        ind2=ts<thre2(j);
        mask1=img2dto3dmask([a1,a2,a3],dim(ind1,:));
        mask2=img2dto3dmask([a1,a2,a3],dim(ind2,:));
        cc1=bwconncomp(mask1,18);
        cc2=bwconncomp(mask2,18);
        ss1=max(cellfun(@length,cc1.PixelIdxList));
        ss2=max(cellfun(@length,cc2.PixelIdxList));
        if isempty(ss1)
            maxcl1(i,j)=0;
        else
            maxcl1(i,j)=ss1;
        end
        if isempty(ss2)
            maxcl2(i,j)=0;
        else
            maxcl2(i,j)=ss2;
        end
    end
    i
end

%real cluster
ts=BWAS_Tregression(design,X);%t value
positive_cluster=cell([mm,1]);  
negative_cluster=cell([mm,1]);
for j=1:mm     %%% j is number of CDT
    positive_cluster{j,2}=thre1(j);
    negative_cluster{j,2}=thre2(j);
    positive_cluster{j,3}=CDT(j);
    negative_cluster{j,3}=CDT(j);
    ind1=ts>thre1(j);
    ind2=ts<thre2(j);
    mask1=img2dto3dmask([a1,a2,a3],dim(ind1,:));
    mask2=img2dto3dmask([a1,a2,a3],dim(ind2,:));
    cc1=bwconncomp(mask1,18);
    cc2=bwconncomp(mask2,18);
    ss1=sort(cellfun(@length,cc1.PixelIdxList),'descend');%Sorts in descending order
    ss2=sort(cellfun(@length,cc2.PixelIdxList),'descend');
    if ~isempty(ss1)
        for jj=1:length(ss1)
            positive_cluster{j,1}(jj,1)=ss1(jj);
            positive_cluster{j,1}(jj,2)=mean(maxcl1(:,j)>ss1(jj));
        end
    else
        positive_cluster{j,1}=[];
    end
    if ~isempty(ss2)
        for jj=1:length(ss2)
            negative_cluster{j,1}(jj,1)=ss2(jj);
            negative_cluster{j,1}(jj,2)=mean(maxcl2(:,j)>ss2(jj));
        end
    else
        negative_cluster{j,1}=[];
    end
end

save dyssomnia_perm5000 positive_cluster negative_cluster ts mask1 mask2;

