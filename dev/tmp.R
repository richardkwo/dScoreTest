library(data.table)
library(plyr)

mu.f <- function(x) {
    5 * exp(x)
}

gen_data <- function(n=500, tau=0) {
    X <- matrix(rnorm(n * 3), nrow=n)
    X[,2] <- X[,2] + 0.2 * X[,1]
    X[,3] <- X[,3] + 0.5 * X[,2]
    y <- mu.f(X[,1]) + exp((X[,2]-1)^2 / 4) * tau / sqrt(n) + rnorm(n) * 3
    return(list(X=X, y=y))
}

dat <- gen_data(500, tau=5)
X <- dat$X
y <- dat$y
dat.df <- data.frame(cbind(X,y))
plot(X[,1], y, pch=20, cex=0.4)
points(X[,2], y, pch=20, col="grey", cex=0.4)
points(X[,3], y, pch=20, col="blue", cex=0.4)
curve(mu.f(x), -3, 3, add=TRUE, col="blue", lwd=1)

fit.0 <- glm(y ~ X, family = gaussian(link = "log"), 
             start = c(log(mean(y)), rep(0, 3)))
print(fit.0)

# methods ---------
score_fun <- function(fit, y, X) {
    mu <- predict(fit, newdata=list(X=X), type="response")
    mu * (mu - y)
}

weight_fun <- function(fit, X) {
    mu <- predict(fit, newdata=list(X=X), type="response")
    mu^2
}

predict_fun <- function(fit, X) {
    predict(fit, newdata=list(X=X))
}

fit_method <- function(y, X) {
    glm(y ~ X, family = gaussian(link = "log"), 
        start = c(log(mean(y[y>0])), rep(0, 3)))
}

wls_method <- function(y, X, w) {
    glm(y ~ X, family = gaussian(), weights=w)
}

# test ---- 

dat <- gen_data(500, tau=0)
X <- dat$X
y <- dat$y
fit.0 <- fit_method(y, X)
print(fit.0)

st <- dScoreTest(y, X, score_fun, weight_fun, fit_method, wls_method, 
                 splits = c(1,1,1),
           hunt.style="vanilla", predict_fun=predict_fun, verbose = TRUE)
print(st)
plot(st)

# -------
dat <- gen_data(500, tau=4)
dat.df <- data.frame(y=dat$y, dat$X)
glm.fit <- glm(y ~ ., family = gaussian(link = "log"), 
               data=dat.df, start = rep(1, 4))

# -------
pval.df <- rdply(200, {
    dat <- gen_data(500, tau=0)
    X <- dat$X
    y <- dat$y
    st <- dScoreTest(y, X, score_fun, weight_fun, fit_method, wls_method, 
                       splits = c(1,1),
                       hunt.style="optimal", 
                       predict_fun=predict_fun, 
                       trim.outlier.hunt = TRUE,
                       verbose = FALSE)
    st.1 <- dScoreTest(y, X, score_fun, weight_fun, fit_method, wls_method, 
                       splits = c(1,1),
                       hunt.style="wls", 
                       predict_fun=predict_fun, 
                       trim.outlier.hunt = TRUE,
                       verbose = FALSE)
    st.2 <- dScoreTest(y, X, score_fun, weight_fun, fit_method, wls_method, 
                     splits = c(1,1),
                     hunt.style="vanilla", 
                     predict_fun=predict_fun, 
                     trim.outlier.hunt = TRUE,
                     verbose = FALSE)
    c(opt=st$p.val, wls=st.1$p.val, van=st.2$p.val)
}, .progress = "text")

plot(ecdf(pval.df$opt), cex=0.5, col="orange")
plot(ecdf(pval.df$wls), cex=0.5, add=TRUE, col="blue")
plot(ecdf(pval.df$van), cex=0.5, add=TRUE, col="grey")
abline(a=0, b=1, col="red", lwd=2)
legend("bottomright", legend=c("opt","wls","vanilla"), 
       pch=c(20,20,20), col=c("orange","blue","grey"))


# -------
pval.df <- rdply(100, {
    dat <- gen_data(1000, tau=5)
    dat.df <- data.frame(y=dat$y, dat$X)
    
    glm.fit <- glm(y ~ ., family = gaussian(link = "log"), 
                   data=dat.df, start = rep(1, 4))
    
    st <- gof_test(glm.fit, splits = c(1,1,1))
    st.1 <- gof_test(glm.fit, splits = c(1,1,1), hunt.style = "wls")
    st.2 <- gof_test(glm.fit, splits = c(1,1,1), hunt.style = "vanilla")
    c(opt=st$p.val, wls=st.1$p.val, van=st.2$p.val)
}, .progress = "text")

plot(ecdf(pval.df$opt), cex=0.5, col="orange")
plot(ecdf(pval.df$wls), cex=0.5, add=TRUE, col="blue")
plot(ecdf(pval.df$van), cex=0.5, add=TRUE, col="grey")
abline(a=0, b=1, col="red", lwd=2)
legend("bottomright", legend=c("opt","wls","vanilla"), 
       pch=c(20,20,20), col=c("orange","blue","grey"))
