# Comparison of Modal Regression with Local Regression

source("C:/Users/ARKAJYOTI/Desktop/IITK SEM 4/Robust Statistical Methods/RobustStats/NPmodalReg/R/meanshift1d.R")

## generate data

n = 500
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

par(mfrow=c(2,2))
plot(X, Y.s, ylim = c(-.6, .9))

Xgrid = seq(from = min(X), to = max(X), length.out = 100)
Ygrid = seq(from = min(Y.s), to = max(Y.s), length.out = 5)
mesh = expand.grid(Xgrid, Ygrid)

#modal regression
y.modalreg = modalReg(Xdata = X, Ydata = Y.s, xgrid = mesh[, 1],
                      ygrid = mesh[, 2], hx = h0, hy = h0)

lines(mesh[, 1][1:100], y.modalreg[1:100], lwd = 3, col = 'blue')
lines(mesh[, 1][401:500], y.modalreg[401:500], lwd = 3, col = 'blue')

MS_fit = modalReg(X, Y.s, hx = h0, hy = h0)
#predicted Y on the modal curves of each data point

alpha=0.95

MS_ps = quantile(abs(MS_fit-Y.s), alpha)
# the quantile of loss; this is larger than the Hausdorff distance so that the prediction set has proper coverage

lines(mesh[, 1][1:100], y.modalreg[1:100] + MS_ps, lwd = 3, col = "orange")
lines(mesh[, 1][401:500], y.modalreg[401:500] + MS_ps, lwd=3, col = "orange")
lines(mesh[, 1][1:100], y.modalreg[1:100] - MS_ps, lwd = 3, col = "orange")
lines(mesh[, 1][401:500], y.modalreg[401:500] - MS_ps, lwd = 3, col = "orange")
legend("bottomleft", c("Modal Regression",
                       paste(100*alpha, "% PS, Modal", sep="")),
       lwd = c(5,5), col=c("blue","orange"), cex = 0.6 )

#Comparison to local regression

span_seq = seq(from = 0.1, to = 0.9, length.out = 100)
loc_PS_seq = rep(NA, 100)

for(i in 1:100){
  fit_loc_tmp = loess(Y.s ~ X, span = span_seq[i])
  loc_PS_seq[i] = quantile(fit_loc_tmp$res, alpha)
}
s_opt = span_seq[which(loc_PS_seq == min(loc_PS_seq))]
#pick the optimal span

fit_loc = loess(Y.s~X,span=s_opt)

par(mfrow=c(2,2))

plot(X, Y.s, ylim = c(-.6, .9))
lines(X[order(X)], fit_loc$fitted[order(X)], lwd=5, col="red")
lines(X[order(X)], fit_loc$fitted[order(X)]+quantile(fit_loc$res,alpha), lwd=3, col="orange")
lines(X[order(X)], fit_loc$fitted[order(X)]-quantile(fit_loc$res,alpha), lwd=3, col="orange")
legend("bottomleft", c("Local Regression",paste(100*alpha,"% PS, Local", sep="")),
       lwd=c(5,5), col=c("red","orange") , cex = 0.7)

plot(X, Y.s, ylim = c(-.6, .9))
lines(mesh[, 1][1:100], y.modalreg[1:100], lwd = 3, col = 'blue')
lines(mesh[, 1][401:500], y.modalreg[401:500], lwd = 3, col = 'blue')
lines(mesh[, 1][1:100], y.modalreg[1:100] + MS_ps, lwd = 3, col = "dodgerblue")
lines(mesh[, 1][401:500], y.modalreg[401:500] + MS_ps, lwd=3, col = "dodgerblue")
lines(mesh[, 1][1:100], y.modalreg[1:100] - MS_ps, lwd = 3, col = "dodgerblue")
lines(mesh[, 1][401:500], y.modalreg[401:500] - MS_ps, lwd = 3, col = "dodgerblue")
legend("bottomleft", c("Modal Regression",
                       paste(100*alpha, "% PS, Modal", sep="")),
       lwd = c(5,5), col=c("blue","dodgerblue"), cex = 0.6 )








