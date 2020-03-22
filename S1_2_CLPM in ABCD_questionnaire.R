# %% Step 1.2: CLPM between Sleep and ADHD in ABCD
# P-value was tested by a multi-level block permutation (family relatedness based on questionnaire)
# written by Dr Qiang Luo and modifiedy by Chun Shen
# Email: qluo@fudan.edu.cn
# released on 21 Mar 2020
# please cite: Shen, et al. Biological Psychiatry 2020

require(lavaan)
require(tidyverse)
require(ggplot2)
require(semPlot)
require(ppcor)
require(xlsx)
require(parallel)

# load data
data <- read.xlsx(file="long_data_0815_formplus.xlsx", header=T, sheetIndex = 1)  # need the data file for ABCD

data[data==-999] <- NA

# extract data and rename the columns
data.clpm <- data[,c(1:2, 3,6, 4,7,5,8,  17, 9,13,  10,14,  11,15, 12,16, 18:37, 38:40, 41,44, 42,45, 43,46)]
colnames(data.clpm) <- c("adhd1","adhd2", 
                           "para1","para2",
                           "dyss1","dyss2",
                           "tot1","tot2",
                           "sex",
                           "bmi1","bmi2",
                           "pub1", "pub2",
                           "ses11","ses21",
                           "ses12","ses22",
                           "s1","s2","s3","s4","s5","s6","s7","s8","s9","s10",
                           "s11","s12","s13","s14","s15","s16","s17","s18","s19","s20",
                           "race1", "race2", "race3",
                           "tr12","tr22",
                           "tr13","tr23",
                           "tr14","tr24") # x ~ adhd; y ~ sleep; sex; bmi; pub; ses(dummy 2); site(dummy 20); race; tr(dummy 3);

print(colnames(data)[c(1:2, 3,6, 4,7,5,8,  17, 9,13,  10,14,  11,15, 12,16, 18:37, 38:40, 41,44, 42,45, 43,46)])

# covariates at the baseline: bmi1,pub1,ses11,ses12 sex,s1-s20,race1,race2,race3,tr1_2,tr1_3,tr1_4
# covariates at the followup: bmi2,pub2,ses21,ses22 sex,s1-s20,race1,race2,race3,tr2_2,tr2_3,tr2_4

##########################
# Model 1: cross-lagged path model
##########################


## model 1.1 unconstrained
clpmModel <- 
  '
#Note, the data contain x1-2 and y1-2
#Latent mean Structure with intercepts

kappa =~ 1*x1 + 1*x2
omega =~ 1*y1 + 1*y2

x1 ~ mu1*1 #intercepts
x2 ~ mu2*1
y1 ~ pi1*1
y2 ~ pi2*1


kappa ~~ 0*kappa #variance
omega ~~ 0*omega #variance
kappa ~~ 0*omega #covariance

#laten vars for AR and cross-lagged effects
p1 =~ 1*x1 #each factor loading set to 1
p2 =~ 1*x2
q1 =~ 1*y1
q2 =~ 1*y2


#regressions
p2 ~ alpha2*p1 + beta2*q1 + sex + bmi2 + pub2 + ses21 + ses22 + s1+s2+s3+s4+s5+s6+s7+s8+s9+s10+s11+s12+s13+s14+s15+s16+s17+s18+s19+s20 + race1+race2+race3 + tr22+tr23+tr24 
q2 ~ delta2*q1 + gamma2*p1+ sex + bmi2 + pub2 + ses21 + ses22 + s1+s2+s3+s4+s5+s6+s7+s8+s9+s10+s11+s12+s13+s14+s15+s16+s17+s18+s19+s20 + race1+race2+race3 + tr22+tr23+tr24 
p1 ~ sex + bmi1 + pub1 + ses11 + ses12 + s1+s2+s3+s4+s5+s6+s7+s8+s9+s10+s11+s12+s13+s14+s15+s16+s17+s18+s19+s20 + race1+race2+race3 + tr12+tr13+tr14
q1 ~ sex + bmi1 + pub1 + ses11 + ses12 + s1+s2+s3+s4+s5+s6+s7+s8+s9+s10+s11+s12+s13+s14+s15+s16+s17+s18+s19+s20 + race1+race2+race3 + tr12+tr13+tr14


p1 ~~ p1 #variance
p2 ~~ u2*p2
q1 ~~ q1 #variance
q2 ~~ v2*q2


p1 ~~ q1 #p1 and q1 covariance
p2 ~~ q2'


colnames.original <- colnames(data.clpm)


colnames(data.clpm) <- colnames.original
colnames(data.clpm)[c(1,2,7,8)] <- c("x1", "x2", "y1", "y2")  # x -- adhd; y -- tot
fit.clpmModel.tot <- lavaan(clpmModel, data = data.clpm,
                            missing = 'ML', #for the missing data!
                            int.ov.free = F,
                            int.lv.free = F,
                            auto.fix.first = F,
                            auto.fix.single = F,
                            auto.cov.lv.x = F,
                            auto.cov.y = F,
                            auto.var = F)


colnames(data.clpm) <- colnames.original
colnames(data.clpm)[c(1,2,5,6)] <- c("x1", "x2", "y1", "y2")  # x -- adhd; y -- dyss
fit.clpmModel.dyss <- lavaan(clpmModel, data = data.clpm,
                             missing = 'ML', #for the missing data!
                             int.ov.free = F,
                             int.lv.free = F,
                             auto.fix.first = F,
                             auto.fix.single = F,
                             auto.cov.lv.x = F,
                             auto.cov.y = F,
                             auto.var = F)

colnames(data.clpm) <- colnames.original
colnames(data.clpm)[c(1,2,3,4)] <- c("x1", "x2", "y1", "y2")  # x -- adhd; y -- para
fit.clpmModel.para <- lavaan(clpmModel, data = data.clpm,
                            missing = 'ML', #for the missing data!
                            int.ov.free = F,
                            int.lv.free = F,
                            auto.fix.first = F,
                            auto.fix.single = F,
                            auto.cov.lv.x = F,
                            auto.cov.y = F,
                            auto.var = F)



tot.beta2.z <- 0
tot.gamma2.z <- 0
dyss.beta2.z <- 0
dyss.gamma2.z <- 0
para.beta2.z <- 0
para.gamma2.z <- 0

para.all <- data.frame(tot.beta2.z, tot.gamma2.z, dyss.beta2.z, dyss.gamma2.z, para.beta2.z, para.gamma2.z)

result.tot <- summary(fit.clpmModel.tot, standardized = T, fit.measures = TRUE)
#para.all$tot.beta2 <- result.tot$PE[which(result.tot$PE$label=='beta2'),"std.all"]
#para.all$tot.gamma2 <- result.tot$PE[which(result.tot$PE$label=='gamma2'),"std.all"]
para.all$tot.beta2.z <- result.tot$PE[which(result.tot$PE$label=='beta2'),"z"]
para.all$tot.gamma2.z <- result.tot$PE[which(result.tot$PE$label=='gamma2'),"z"]

result.dyss <- summary(fit.clpmModel.dyss, standardized = T, fit.measures = TRUE)
#para.all$dyss.beta2 <- result.dyss$PE[which(result.dyss$PE$label=='beta2'),"std.all"]
#para.all$dyss.gamma2 <- result.dyss$PE[which(result.dyss$PE$label=='gamma2'),"std.all"]
para.all$dyss.beta2.z <- result.dyss$PE[which(result.dyss$PE$label=='beta2'),"z"]
para.all$dyss.gamma2.z <- result.dyss$PE[which(result.dyss$PE$label=='gamma2'),"z"]

result.para <- summary(fit.clpmModel.para, standardized = T, fit.measures = TRUE)
#para.all$para.beta2 <- result.para$PE[which(result.para$PE$label=='beta2'),"std.all"]
#para.all$para.gamma2 <- result.para$PE[which(result.para$PE$label=='gamma2'),"std.all"]
para.all$para.beta2.z <- result.para$PE[which(result.para$PE$label=='beta2'),"z"]
para.all$para.gamma2.z <- result.para$PE[which(result.para$PE$label=='gamma2'),"z"]


print(para.all)


colnames(data.clpm) <- colnames.original



#p.fdr <- as.data.frame(p.adjust(result$PE[c(40:59),]$pvalue, method  = "BH"))
#rownames(p.fdr) <- result$PE[c(40:59),]$label
#colnames(p.fdr) <- 'p.fdr'
#print(p.fdr)

######
## permutation
######

permID <- read.csv("Pset4crosslagABCD.csv", header = T)
colnames(permID)[1] <- "ID"

cl = makeCluster(12)

# check it works.
clusterEvalQ(cl, runif(10))    # -> this works
clusterSetRNGStream(cl, 100)    # to get reproduciable results

rm(results)
nperm <- 5000
processInput  <- function(i){
  
  # permute
  data.clpm.perm <- data.clpm
  data.clpm.perm[,c(1,2)] <- data.clpm.perm[permID[,i+1],c(1,2)]  # permute ADHD score at the baseline
  # model for dyss
  colnames(data.clpm.perm) <- colnames.original
  colnames(data.clpm.perm)[c(1,2,5,6)] <- c("x1", "x2", "y1", "y2")  # x -- adhd; y -- dyss
  fit.clpmModel.dyss.perm <- lavaan(clpmModel, data = data.clpm.perm,
                                    missing = 'ML', #for the missing data!
                                    int.ov.free = F,
                                    int.lv.free = F,
                                    auto.fix.first = F,
                                    auto.fix.single = F,
                                    auto.cov.lv.x = F,
                                    auto.cov.y = F,
                                    auto.var = F)
  result.dyss <- summary(fit.clpmModel.dyss.perm, standardized = T, fit.measures = TRUE)
  dyss.beta2.z <- result.dyss$PE[which(result.dyss$PE$label=='beta2'),"z"]
  dyss.gamma2.z <- result.dyss$PE[which(result.dyss$PE$label=='gamma2'),"z"]
  
  
  
  # model for para
  colnames(data.clpm.perm) <- colnames.original
  colnames(data.clpm.perm)[c(1,2,3,4)] <- c("x1", "x2", "y1", "y2")  # x -- adhd; y -- para
  fit.clpmModel.para.perm <- lavaan(clpmModel, data = data.clpm.perm,
                                    missing = 'ML', #for the missing data!
                                    int.ov.free = F,
                                    int.lv.free = F,
                                    auto.fix.first = F,
                                    auto.fix.single = F,
                                    auto.cov.lv.x = F,
                                    auto.cov.y = F,
                                    auto.var = F)
  result.para <- summary(fit.clpmModel.para.perm, standardized = T, fit.measures = TRUE)
  para.beta2.z <- result.para$PE[which(result.para$PE$label=='beta2'),"z"]
  para.gamma2.z <- result.para$PE[which(result.para$PE$label=='gamma2'),"z"]
  
  
  
  # model for tot
  colnames(data.clpm.perm) <- colnames.original
  colnames(data.clpm.perm)[c(1,2,7,8)] <- c("x1", "x2", "y1", "y2")  # x -- adhd; y -- tot
  fit.clpmModel.tot.perm <- lavaan(clpmModel, data = data.clpm.perm,
                                   missing = 'ML', #for the missing data!
                                   int.ov.free = F,
                                   int.lv.free = F,
                                   auto.fix.first = F,
                                   auto.fix.single = F,
                                   auto.cov.lv.x = F,
                                   auto.cov.y = F,
                                   auto.var = F)
  result.tot <- summary(fit.clpmModel.tot.perm, standardized = T, fit.measures = TRUE)
  tot.beta2.z <- result.tot$PE[which(result.tot$PE$label=='beta2'),"z"]
  tot.gamma2.z <- result.tot$PE[which(result.tot$PE$label=='gamma2'),"z"]
  
  
  para.perm.all <- data.frame(tot.beta2.z, tot.gamma2.z, dyss.beta2.z, dyss.gamma2.z, para.beta2.z, para.gamma2.z)
  
  return(para.perm.all)
}
clusterEvalQ(cl, library(lavaan))  # load the library to every core
clusterExport(cl, c("processInput","data.clpm", "permID", "colnames.original", "clpmModel")) # export the variables to every core
results <- parLapply(cl, 1:nperm, function(i) processInput(i))
para.perm.all <- data.frame(matrix(0,nrow = nperm, ncol = 6))
colnames(para.perm.all) <- colnames(para.all) 
for (i in c(1:nperm)){
  para.perm.all[i,] <- as.data.frame(results[[i]])
}
rm(results)
max.tot  <- matrix(apply(abs(para.perm.all[,c(1,2)]), 1, function(x){ max(x) }), nrow = nperm, ncol = 1)
max.dimension <- matrix(apply(abs(para.perm.all[,c(3,6)]), 1, function(x){ max(x) }), nrow = nperm, ncol = 1)
perm.p <- matrix(1, nrow = 1, ncol = 6)
perm.p[1,c(1,2)] <- apply(abs(para.all[c(1:2)]), 2, function(x) {sum(max.tot > x) / nperm})
perm.p[1,c(3:6)] <- apply(abs(para.all[c(3:6)]), 2, function(x) {sum(max.dimension > x) / nperm})
colnames(perm.p) <- colnames(para.all)
print(perm.p)


stopCluster(cl)
