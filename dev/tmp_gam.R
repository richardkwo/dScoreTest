library(mgcv)
set.seed(42)
dat <- gamSim(eg=1, n=400, dist="normal", scale=2, verbose = FALSE)
dat.0 <- dat[,1:5]
# well-specified
fit.0 <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),data=dat.0)
test.0 <- gof_test(fit.0)
# also well-specified
fit.1 <- gam(y~s(x0)+s(x1)+s(x2),data=dat.0)
test.1 <- gof_test(fit.1)
plot(test.1)
# misspecified
dat.1 <- dat.0
dat.1$y <- dat.1$y * dat$f0 
fit.2 <- gam(y~s(x0)+s(x1)+s(x2)+s(x3), data=dat.1)
test.2 <- gof_test(fit.2)
plot(test.2)

# wls -------

dat <- gamSim(eg=1, n=400, dist="normal", scale=2, verbose = FALSE)
dat.0 <- dat[,1:5]
# well-specified
fit.0 <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),data=dat.0)
gof_test(fit.0, hunt.style = "wls")
# also well-specified
fit.1 <- gam(y~s(x0)+s(x1)+s(x2),data=dat.0)
gof_test(fit.1, hunt.style = "wls")
# misspecified
dat.1 <- dat.0
dat.1$y <- dat.1$y * dat$f0^2 
fit.2 <- gam(y~s(x0)+s(x1)+s(x2)+s(x3), data=dat.1)
test.2 <- gof_test(fit.2, hunt.style = "wls")
plot(test.2)



