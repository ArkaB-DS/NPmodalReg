# Based on Chapter 7: Density Estimation: Erupting Geysers and Star Clusters


## kernel functions
rec <- function(x) (abs(x) < 1) * 0.5
tri <- function(x) (abs(x) < 1) * (1 - abs(x))

x <- seq(from = -3, to = 3, by = 0.001)

## Figure 7.1 Three commonly used kernel functions.
plot(x, rec(x), type = "l", ylim = c(0, 1),
     ylab = expression(K(x)), col = "seagreen",
     main = "Three commonly used kernel functions")
lines(x, tri(x), lty = 2, col = "blue")
lines(x, dnorm(x), lty = 3, col = "red")
legend("topright", legend = c("Rectangular", "Triangular", "Gaussian"),
       lty = 1:3, col = c("seagreen", "blue", "red"))


## artificial data points
x <- c(0, 1, 1.1, 1.5, 1.9, 2.8, 2.9, 3.5)
n <- length(x)

## create a grid
xgrid <- seq(from = min(x) - 1, to = max(x) + 1, by = 0.01)

## window width
h <- 0.4

# create bumps with width h and shape "Gaussian"
bumps <- sapply(x, function(k) dnorm((xgrid - k)/h)/(n * h))

## Fig 7.2: Kernel estimate showing the contributions of Gaussian kernels evaluated
## for the individual observations with bandwidth h = 0.4.
plot(xgrid, rowSums(bumps), ylab = expression(hat(f)(x)),
        type = "l", xlab = "x", lwd = 2, col = "red", 
     main = "Kernel estimate showing the contributions of Gaussian kernels evaluated
for the individual observations with bandwidth h = 0.4.")
rug(x, lwd = 2)
out <- apply(bumps, 2, function(b) lines(xgrid, b, col = "blue"))

## Epanechnikov kernel 
epa <- function(x, y) ((x^2 + y^2) < 1) * 2/pi * (1 - x^2 - y^2)

x <- seq(from = -1.1, to = 1.1, by = 0.05)

## Fig 7.3 Epanechnikov kernel for a grid between (−1.1,−1.1) and (1.1, 1.1).
epavals <- sapply(x, function(a) epa(a, x))
persp(x = x, y = x, z = epavals, xlab = "x", ylab = "y",
      zlab = expression(K(x, y)), theta = -35, phi = 15,
      axes = TRUE, box = TRUE, 
      main = "Epanechnikov kernel for a grid between (−1.1,−1.1) and (1.1, 1.1).",
      col = "lightblue")

## Density estimates of the geyser eruption data imposed on a histogram
## of the data.
x <- faithful$waiting

## Figure 7.4 
par(mfrow = c(1, 3))
hist(x, xlab = "Waiting times (in min.)", ylab = "Relative Frequency",
     probability = TRUE, main = "Gaussian kernel",
     border = "gray")
lines(density(x, width = 12), lwd = 2)
rug(x)
hist(x, xlab = "Waiting times (in min.)", ylab = "Relative Frequency",
     probability = TRUE, main = "Rectangular kernel",
     border = "gray")
lines(density(x, width = 12, window = "rectangular"), lwd = 2)
rug(x)
hist(x, xlab = "Waiting times (in min.)", ylab = "Relative Frequency",
     probability = TRUE, main = "Triangular kernel",
     border = "gray")
lines(density(x, width = 12, window = "triangular"), lwd = 2)
rug(x)


## 7.3.1 A Parametric Density Estimate for the Old Faithful Data

logL <- function(param, x) {
  d1 <- dnorm(x, mean = param[2], sd = param[3])
  d2 <- dnorm(x, mean = param[4], sd = param[5])
  -sum(log(param[1] * d1 + (1 - param[1]) * d2))
  }
startparam <- c(p = 0.5, mu1 = 50, sd1 = 3, mu2 = 80, sd2 = 3)
opp <- optim(startparam, logL, x = faithful$waiting,
                method = "L-BFGS-B",
                lower = c(0.01, rep(1, 4)),
                upper = c(0.99, rep(200, 4)))
opp
# $par
# p        mu1        sd1        mu2        sd2 
# 0.3608912 54.6121380  5.8723776 80.0934106  5.8672829 
# 
# $value
# [1] 1034.002
# 
# $counts
# function gradient 
# 55       55 
# 
# $convergence
# [1] 0
# 
# $message
# [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"

opar <- as.list(opp$par)
rx <- seq(from = 40, to = 110, by = 0.1)
d1 <- dnorm(rx, mean = opar$mu1, sd = opar$sd1)
d2 <- dnorm(rx, mean = opar$mu2, sd = opar$sd2)
f <- opar$p * d1 + (1 - opar$p) * d2

## Figure 7.7 Fitted normal density and two-component normal mixture for geyser
## eruption data.
hist(x, probability = TRUE, xlab = "Waiting times (in min.)",
        border = "gray", xlim = range(rx), ylim = c(0, 0.06),
        main = "")
lines(rx, f, lwd = 2)
lines(rx, dnorm(rx, mean = mean(x), sd = sd(x)), lty = 2,
         lwd = 2)
legend(50, 0.06, lty = 1:2, bty = "n",
          legend = c("Fitted two-component mixture density",
                       "Fitted single normal density"))

library("mclust")
mc <- Mclust(faithful$waiting)
mc
mc$parameters$mean
sqrt(mc$parameters$variance$sigmasq)

library("flexmix")
fl <- flexmix(waiting ~ 1, data = faithful, k = 2)
parameters(fl, component = 1)
parameters(fl, component = 2)

library("boot")
fit <- function(x, indx) {
  a <- Mclust(x[indx], minG = 2, maxG = 2)$parameters
  if (a$pro[1] < 0.5)
    return(c(p = a$pro[1], mu1 = a$mean[1],
             mu2 = a$mean[2]))
    return(c(p = 1 - a$pro[1], mu1 = a$mean[2],
             mu2 = a$mean[1]))
}

bootpara <- boot(faithful$waiting, fit, R = 1000)
boot.ci(bootpara, type = "bca", index = 1)
boot.ci(bootpara, type = "bca", index = 2)
boot.ci(bootpara, type = "bca", index = 3)
bootplot <- function(b, index, main = "") {
  dens <- density(b$t[,index])
  ci <- boot.ci(b, type = "bca", index = index)$bca[4:5]
  est <- b$t0[index]
  plot(dens, main = main)
  y <- max(dens$y) / 10
  segments(ci[1], y, ci[2], y, lty = 2)
  points(ci[1], y, pch = "(")
  points(ci[2], y, pch = ")")
  points(est, y, pch = 19)
}
## Figure 7.8 Bootstrap distribution and confidence intervals for the mean estimates
## of a two-component mixture for the geyser data.
layout(matrix(1:2, ncol = 2))
bootplot(bootpara, 2, main = expression(mu[1]))
bootplot(bootpara, 3, main = expression(mu[2]))


## A contour plot of the bivariate density estimate of the CYGOB1 data, i.e., 
## a two-dimensional graphical display for a three-dimensional problem.

library("KernSmooth")

## Energy output and surface termperature for Star Cluster CYG OB1
data("CYGOB1", package = "HSAUR")

## lookup ??dpik, ??bkde2D 
CYGOB1d <- bkde2D(CYGOB1, bandwidth = sapply(CYGOB1, dpik))
par(mfrow = c(1, 2))
## Fig 7.5
contour(x = CYGOB1d$x1, y = CYGOB1d$x2, z = CYGOB1d$fhat,
           xlab = "log surface temperature",
           ylab = "log light intensity")
image(x = CYGOB1d$x1, y = CYGOB1d$x2, z = CYGOB1d$fhat,
        xlab = "log surface temperature",
        ylab = "log light intensity")

## Figure 7.6 The bivariate density estimate of the CYGOB1 data, here shown in a
## three-dimensional fashion using the persp function.
par(bg = "seagreen")
persp(x = CYGOB1d$x1, y = CYGOB1d$x2, z = CYGOB1d$fhat,
      xlab = "log surface temperature",
      ylab = "log light intensity",
      zlab = "estimated density",
      theta = -35, phi = 10, axes = TRUE, box = TRUE,
      col = "magenta", bg = "black")
