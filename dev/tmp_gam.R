library(mgcv)
library(dScoreTest)

set.seed(42)
# prostate cancer data
prostate.dat <- read.table("https://hastie.su.domains/ElemStatLearn/datasets/prostate.data", 
                       header = TRUE)

prostate.dat$train <- NULL
colnames(prostate.dat)[c(1,2,4,6,9)] <- c("log.can.vol", 
                                          "log.weight", 
                                          "log.BPH", 
                                          "log.cap.pen", 
                                          "log.PSA")

# lpsa (log PSA) is the outcome
hist(prostate.dat$log.PSA)

# svi (seminal vesicle invasion) is binary
table(prostate.dat$svi)

# Gleason score is categorical
table(prostate.dat$gleason)

# so only put linear terms for svi and gleason later

# GAM ------

# Model and visualize how log.PSA varies according to log(can.vol) and log(weight) -----

# fit a model without interaction: 
fit.0 <- gam(log.PSA ~ ti(log.can.vol, k=6) + 
                 ti(log.weight, k=6) + 
                 s(x=age,k=7,fx=F,bs="cr",m=2)+
                 s(x=log.BPH,k=7,fx=F,bs="cr",m=2) + 
                 s(x=log.cap.pen,k=7,fx=F,bs="cr",m=2)+
                 s(x=pgg45,k=7,fx=F,bs="cr",m=2) + 
                 svi + gleason, 
             data=prostate.dat, method="GCV.Cp")

# a model with interaction
fit.1 <- gam(log.PSA ~ ti(log.can.vol, k=6) + 
                 ti(log.weight, k=6) + 
                 ti(log.can.vol, log.weight, k=c(6,6)) + 
                 s(x=age,k=7,fx=F,bs="cr",m=2)+
                 s(x=log.BPH,k=7,fx=F,bs="cr",m=2) + 
                 s(x=log.cap.pen,k=7,fx=F,bs="cr",m=2)+
                 s(x=pgg45,k=7,fx=F,bs="cr",m=2) + 
                 svi + gleason, 
             data=prostate.dat, method="GCV.Cp")

# visualize
split.screen(c(1,2))
screen(1)
vis.gam(fit.0, view=c("log.weight", "log.can.vol"), theta=35,phi=25, 
        zlab="fitted log.PSA", main="without interaction")
screen(2)
vis.gam(fit.1, view=c("log.weight", "log.can.vol"), theta=35,phi=25, 
        zlab="fitted log.PSA", main="with interaction")
close.screen(all = TRUE)

# Is the with-interaction GAM model is well-specified?

gof.1 <- gof_test(fit.1)
gof.1

plot(gof.1)

# no indication of misspecification

# now test with or without interaction ---
compare_models(fit.0, fit.1)

# compare this with mgcv's approximate significance test for the interaction
summary(fit.1)
# which gives p-val 0.015 for the interaction

# sample size is small --- p.val can be sensitive to data split
replicate(10, compare_models(fit.0, fit.1)$p.val)

# run 100 times and take a heavy-tailed combination (e.g., harmonic mean)
pvals.exch <- replicate(100, compare_models(fit.0, fit.1)$p.val)
pval.interaction <- 1 / mean(1 / pvals.exch)
