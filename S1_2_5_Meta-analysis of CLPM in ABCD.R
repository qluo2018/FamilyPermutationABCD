# %% Step 1.4: meta-analysis for CLPM in ABCD
# %% written by Shen Chun, cshen17@fudan.edu.cn
# %% reviewed by Dr Qiang Luo, qluo@fudan.edu.cn
# %% released on 21 Mar 2020
# %% please cite: Shen, et al. Biological Psychiatry 2020

library(metafor)

# load standard coefficients by site
total_sleep <- read_excel("CLPM_B_by_site.xlsx")

#ADHD->total sleep
dat <- escalc(measure = "COR", ri = b1, ni=n, data = total_sleep, append = TRUE)
res <- rma(yi, vi, data = dat)
forest(res)
summary(res)


#total sleep->ADHD
dat2 <- escalc(measure = "COR", ri = b2, ni=n, data = total_sleep, append = TRUE)
res2 <- rma(yi, vi, data = dat2)
forest(res2)
summary(res2)

#ADHD->dysomnia
dat3 <- escalc(measure = "COR", ri = b3, ni=n, data = total_sleep, append = TRUE)
#dat3 <- escalc(measure = "UCOR", ri = b3, ni=n, data = total_sleep, vtype = "US",append = TRUE)
res3 <- rma(yi, vi, data = dat3)
forest(res3)
summary(res3)

#dysomnia->ADHD
dat4 <- escalc(measure = "COR", ri = b4, ni=n, data = total_sleep, append = TRUE)
#dat4 <- escalc(measure = "COR", ri = b4, ni=n, data = total_sleep, vtype = "LS", append = TRUE)
res4 <- rma(yi, vi, data = dat4)
forest(res4)
summary(res4)

#ADHD->parasomnia
dat5 <- escalc(measure = "COR", ri = b5, ni=n, data = total_sleep, append = TRUE)
res5 <- rma(yi, vi, data = dat5)
forest(res5)
summary(res5)

#parasomnia->ADHD
dat6 <- escalc(measure = "COR", ri = b6, ni=n, data = total_sleep, append = TRUE)
res6 <- rma(yi, vi, data = dat6)
forest(res6)
summary(res6)
