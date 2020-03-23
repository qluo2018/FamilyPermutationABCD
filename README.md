# FamilyPermutationABCD

Key words:
the ABCD cohort; family relatedness; multi-level block permutation; 
whole-brain voxel-wise mediation analysis; random intercept cross-lagged panel model (RI-CLPM).

Data and Code of paper: 

Shen et al. What is the Link between Attention-Deficit/Hyperactivity Disorder and Sleep Disturbance? 
A multimodal examination of longitudinal relationships and brain structure using large-scale population-based cohorts. 
Biological Psychiatry, 2020. Acceptted

NOTE
1. Some codes were adapted from others. Clear copyright and reference was stated in the relevant scripts.
2. Per the data accessing policy of the GRIP (http://www.gripinfo.ca/Grip/Public/www/QuoiNeuf/en/faitssaillants.asp?langue=en),  
the QLSCD data could be used only after approval. Therefore, the QLSCD data in this folder are not allowed to be shared to anyone else.
3. Per the data accessing policy of the NDA agreement, the usage of the ABCD data needs the approval by NIMH. The ABCD cohort Release v1.0 (DOI: 10.15154/1460410) was used in this study.  
4. The GMV file was preprocessed by a PhD student in our group, Mr Jingnan Du using the pipeline descripted in the Method section of the main text. 
5. If the data and codes are used in your work, please cite the above reference, namely Shen et al. 2020. 

SUMMARY （Matlab2018b; MPlus; R3.5.1; PALM:https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM; KING:http://people.virginia.edu/~wc9c/KING/）

Step 1: Cross-lagged Analysis between ADHD symptoms and sleep disturbances.
   
     Step 1.1: the RI-CLPM analysis of the QLSCD sample; (S1_1_RI-CLPM in QLSCD.inp);
   
     Step 1.2: the CLPM analysis of the ABCD sample:
     
          Step 1.2.1: define the multi-level blocks according to the self-reported family relatedness provided by a questionnaire
          (acspsw02.txt) in the ABCD data release v1, and generate the block permuted sample by calling the PALM function
          (S1_2_1_EBforCLPM_questionnaire.m);
     
          Step 1.2.2: define the multi-level blocks by the genetic kinship estimated using the KING software and generate the permuted
          sample ( S1_2_2_EBforCLPM_genetic.m);
     
          Step 1.2.3 CLPM by the questionnaire-defined multi-level block permutation  (S1_2_3_CLPM in ABCD_questionnaire.R)
     
          Step 1.2.4 CLPM by the genetic-defined multi-level block permutation  (S1_2_4_CLPM in ABCD_gene.R)
     
          Step 1.2.5 meta-analysis of the cross-lagged path coefficients among the data collection sites in ABCD (S1_2_5_Meta-analysis 
          of CLPM in ABCD.R)

Step 2: Whole-brain voxel-wise association analysis between GMV and sleep/ADHD with multi-level block permutation  
(S2_MRI_analysis_with_block_permutation.m)

Step 3: Mediation analysis of GMV-->ADHD-->Sleep; Analysis of ROI; and Exporatory analysis of the whole brain voxels by using the 3M toolboxl (S3_Mediation_overlapping_clusters.m).

Step 4: Transcriptomic analysis
 
    Step 4.1 Preprocessing of the Allen Brain Gene expression (S4_1_AHBA_datapreprocess.m);
  
    Step 4.2 Extracting the t-stats of the mediation for each coordinate as an averaged value within a small volume centered at each
    tissue location (S4_2_mediation_t_value_extract.m);
  
    Step 4.3 Identifying gene associations by the partial least squared regression  (S4_3_PLSregression.m).


Contact: 

    Qiang Luo, PhD, Visiting Fellow at the Clare Hall, Cambridge
    Associate Principal Investigator
    Institute of Science and Technology for Brain-Inspired Intelligence (ISTBI)
    Fudan University
    Email: mrqiangluo@gmail.com;      qluo@fudan.edu.cn
    Office Phone: 86-21-65648454
    https://sites.google.com/site/qluochina/
    http://homepage.fudan.edu.cn/qiangluo/

Team Members:  

    Chun Shen, MSc,  cshen17@fudan.edu.cn
    Jingnan Du, BA,  jingnandu3@gmail.com
    Xingzhong Zhao, BA, 18210850006@fudan.edu.cn

The current version was relesed on 22 Mar 2020. 

Citation of this package:
Shen et al. What is the Link between Attention-Deficit/Hyperactivity Disorder and Sleep Disturbance? 
A multimodal examination of longitudinal relationships and brain structure using large-scale population-based cohorts. 
Biological Psychiatry, 2020. Acceptted

 
