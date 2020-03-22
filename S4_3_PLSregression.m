% PLS analysis, code adapted from Whitaker et al (2016) and Vertes et al.(2016)
% Reference: 
% Whitaker KJ et al. (2016): Adolescence is associated with genomically patterned 
% consolidation of the hubs of the human brain connectome. Proc Natl Acad Sci U S A 113:9105-9110.
% Vertes PE et al. (2016): Gene transcription profiles associated with inter-modular 
% hubs and connection distance in human functional magnetic resonance imaging networks. Philos T R Soc B 371.
% Morgan SE et al (2019): Cortical patterning of abnormal morphometric similarity 
% in psychosis is associated with brain expression of schizophrenia-related
% genes. Proc Natl Acad Sci U S A 116 (19) 9604-9609.
% https://github.com/SarahMorgan/Morphometric_Similarity_SZ

%% written by Shen Chun, cshen17@fudan.edu.cn
%% reviewed by Dr Qiang Luo, qluo@fudan.edu.cn
%% released on 21 Mar 2020
%% please cite: Shen, et al. Biological Psychiatry 2020

% load data (using subcortical region as an example)
load('ab_tvalue_80ROI.mat','tvalue_sub_l'); 
% description of data file:
% tvalue_sub_l ----  the t-stats of the mediation effect given by 3M
% toolbox at each tissue location (averaged in a 6-mm samll ball centered 
% at this location) in the subcortical areas in the Allen Brain 
% tvalue_cort_l ----  the t-stats of the mediation effect given by 3M
% toolbox at each tissue location (averaged in a 6-mm samll ball centered 
% at this location) in the cortical areas in the Allen Brain 
load('AHBA_Mean_scaled_reanote.mat');  
% description of the data
% cort_expMeanScaled   ---- all gene expression data of the cortical tissues 
%                           from 6 donors 
% cort_l_expMS         ---- 784 cortical tissues in the left hemisphere for 
%                           15408 genes
% geneSymbol           ---- gene names  
% genes                ---- entrez gene id
% sub_expMeanScaled    ---- all gene expression data of the subcortical tissues 
%                           from 6 donors 
% sub_l_expMS          ---- 182 subcortical tissues in the left hemisphere for 
%                           15408 genes
%%
X=sub_l_expMS; % Predictors
Y=tvalue_sub_l; % Response variable

X=zscore(X);
Y=zscore(Y);

%perform full PLS and plot variance in Y explained by top 15 components
%typically top 2 or 3 components will explain a large part of the variance
%(hopefully!)
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,15);
%PCTVAR containing the percentage of variance explained by the model. 
%The first row of PCTVAR contains the percentage of variance explained in X by each PLS component, 
%and the second row contains the percentage of variance explained in Y
dim=15;
plot(1:dim,cumsum(100*PCTVAR(2,1:dim)),'-o','LineWidth',1.5,'Color',[140/255,0,0]);
set(gca,'Fontsize',14)
xlabel('Number of PLS components','FontSize',14);
ylabel('Percent Variance Explained in Y','FontSize',14);
grid on

dim=15;
plot(1:dim,cumsum(100*PCTVAR(1,1:dim)),'-o','LineWidth',1.5,'Color',[140/255,0,0]);
set(gca,'Fontsize',14)
xlabel('Number of PLS components','FontSize',14);
ylabel('Percent Variance Explained in X','FontSize',14);
grid on

%%% plot correlation of PLS component 1 with t-statistic:
figure
plot(XS(:,1),Y,'r.')% the predictor scores XS, that is, the PLS components that are linear combinations of the variables in X
[R,P]=corrcoef(XS(:,1),Y);
xlabel('XS scores for PLS component 1','FontSize',14);
ylabel('a*b t-statistic','FontSize',14);
grid on
l = lsline;
l.Color = 'k';

% permutation testing to assess significance of PLS result as a function of
% the number of components (dim) included:
rep=5000;
dim = 1;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
Rsquared_real_y = cumsum(100*PCTVAR(2,1:dim));
Rsquared_real_x = cumsum(100*PCTVAR(1,1:dim));

Pvalue_perm_y = ones(1,8);
Pvalue_perm_x = ones(1,8);
for i = 1:dim
    i
    parfor j = 1:rep
    order=randperm(size(Y,1));
    Yp=Y(order,:);
    
    [XL_p,YL_p,XS_p,YS_p,BETA_p,PCTVAR_p,MSE_p,stats_p]=plsregress(X,Yp,dim);
    temp_y=cumsum(100*PCTVAR_p(2,1:dim));
    Rsq_y(j) = temp_y(i);
    temp_x=cumsum(100*PCTVAR_p(1,1:dim));
    Rsq_x(j) = temp_x(i);
    end
    Pvalue_perm_y(i) = length(find(Rsq_y>=Rsquared_real_y(i)))/rep;
    Pvalue_perm_x(i) = length(find(Rsq_x>=Rsquared_real_x(i)))/rep;
end
%%
%Bootstrap to get the gene list:
%in order to estimate the error in estimating each gene's PLS1 weight, then
%to caculate the z score
%number of bootstrap iterations:
bootnum=5000;

% Do PLS in 2 dimensions (with 2 components):
dim=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);

%store regions' IDs and weights in descending order of weight for both components:
[R1,p1]=corr([XS(:,1),XS(:,2)],Y);%XS: predicted X score

%align PLS components with desired direction for interpretability 
if R1(1,1)<0  %this is specific to the data shape we were using - will need ammending
    stats.W(:,1)=-1*stats.W(:,1);
    XS(:,1)=-1*XS(:,1);
end
if R1(2,1)<0 %this is specific to the data shape we were using - will need ammending
    stats.W(:,2)=-1*stats.W(:,2);
    XS(:,2)=-1*XS(:,2);
end

%real weight and sorted by weight
[PLS1w,x1] = sort(stats.W(:,1),'descend');%W: A p-by-ncomp matrix of PLS weights so that XS = X0*W.
PLS1ids=genes(x1);
[PLS2w,x2] = sort(stats.W(:,2),'descend');
PLS2ids=genes(x2);

%define variables for storing the (ordered) weights from all bootstrap runs
PLS1weights=[];
PLS2weights=[];

%start bootstrap
parfor i=1:bootnum
    myresample = randsample(size(X,1),size(X,1),1);%replacement
    res(i,:)=myresample; %store resampling out of interest
    Xr=X(myresample,:); % define X for resampled subjects
    Yr=Y(myresample,:); % define X for resampled subjects
    [XL_b,YL_b,XS_b,YS_b,BETA_b,PCTVAR_b,MSE_b,stats_b]=plsregress(Xr,Yr,dim); %perform PLS for resampled data
      
    temp1=stats_b.W(:,1);%extract PLS1 weights
    newW1=temp1(x1); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS1w,newW1)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW1=-1*newW1;
    end
    PLS1weights=[PLS1weights,newW1];%store (ordered) weights from this bootstrap run
    
    temp2=stats_b.W(:,2);%extract PLS2 weights
    newW2=temp2(x2); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS2w,newW2)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW2=-1*newW2;
    end
    PLS2weights=[PLS2weights,newW2]; %store (ordered) weights from this bootstrap run    
end

%get standard deviation of weights from bootstrap runs
PLS1sw=std(PLS1weights');
PLS2sw=std(PLS2weights');

%get bootstrap weights (Z)
PLS1Z=PLS1w./PLS1sw';
PLS2Z=PLS2w./PLS2sw';

save PLS_sub_Mean_Scaled_Reanote PLS1ids PLS2ids PLS1Z PLS2Z;
