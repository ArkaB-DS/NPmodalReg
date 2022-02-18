source("./R/meanshift1d.R")
library(pracma)

## generate data

n = 200
sigma = 0.25

X1 = runif(n)
X2 = runif(n)
X = c(X1, X2)

Y1 = rnorm(n, 1.5 + 0.5 * sin(X1 * 3 * pi), sd = sigma)
Y2 = rnorm(n, 0.5 * sin(X2 * 3 * pi), sd = sigma)
Y = c(Y1, Y2)

#scaled Y
Y.s = Y/(sd(Y)/sd(X))

#smoothing parameter (non-optimized)
h0 = 0.075

Xgrid = seq(from = min(X), to = max(X), length.out = 50)
Ygrid = seq(from = min(Y.s), to = max(Y.s), length.out = 5)
mesh = expand.grid(Xgrid, Ygrid)

#modal regression
y.modalreg = modalReg(Xdata = X, Ydata = Y.s, xgrid = mesh[, 1],
                      ygrid = mesh[, 2], hx = h0, hy = h0)

alpha=0.95
# confidence set
times = 30
pointwise_CI=numeric(length(mesh[,1]))
for(i in 1:length(mesh[,1]))
{
  cat(i,"\n")
  y.modalregx = modalReg(Xdata = X, Ydata = Y.s, xgrid = rep(mesh[i,1],length(mesh[,1])),
                        ygrid = mesh[, 2], hx = h0, hy = h0)
  haus=numeric(times)
  for(t in 1:times)
  {
    x_boot<-sample(X, length(X), replace = TRUE)
    y_boot<-sample(y.modalreg, length(y.modalreg), replace = TRUE)
    ystar.modalregx = modalReg(Xdata = x_boot, Ydata = y_boot, xgrid = rep(mesh[i,1],length(mesh[,1])),
                           ygrid = mesh[, 2], hx = h0, hy = h0)
    haus[t] <- hausdorff_dist(ystar.modalregx,y.modalregx)
  }
  pointwise_CI[i] = quantile(haus, alpha)
}

uniform_CI=max(pointwise_CI)

plot(X, Y.s, cex =0.5,  ylim = c(-0.5, 1), ylab="Y")
lines(mesh[, 1][1:50], y.modalreg[1:50], lwd = 2, col = 'blue')
lines(mesh[, 1][1:50], y.modalreg[1:50]+uniform_CI, lwd = 2, col = 'red')
lines(mesh[, 1][1:50], y.modalreg[1:50]-uniform_CI, lwd = 2, col = 'red')
lines(mesh[, 1][201:250], y.modalreg[201:250], lwd = 2, col = 'blue',)
lines(mesh[, 1][201:250], y.modalreg[201:250]+uniform_CI, lwd = 2, col = 'red')
lines(mesh[, 1][201:250], y.modalreg[201:250]-uniform_CI, lwd = 2, col = 'red')

plot(X, Y.s, cex =0.5,  ylim = c(-0.5, 1), ylab="Y")
lines(mesh[, 1][1:50], y.modalreg[1:50], lwd = 2, col = 'blue')
lines(mesh[, 1][1:50], y.modalreg[1:50]+pointwise_CI[1:20], lwd = 2, col = 'red')
lines(mesh[, 1][1:50], y.modalreg[1:50]-pointwise_CI[1:20], lwd = 2, col = 'red')
lines(mesh[, 1][201:250], y.modalreg[201:250], lwd = 2, col = 'blue',)
lines(mesh[, 1][201:250], y.modalreg[201:250]+pointwise_CI[81:100], lwd = 2, col = 'red')
lines(mesh[, 1][201:250], y.modalreg[201:250]-pointwise_CI[81:100], lwd = 2, col = 'red')

