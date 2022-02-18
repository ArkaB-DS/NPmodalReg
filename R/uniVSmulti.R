source("./R/meanshift1d.R")
source("./R/unimode.R")

## generate data

n = 500
sigma = 0.25

X1 = runif(n/2)
X2 = runif(n)
X3 = runif(n/2)
X = c(X1, X2, X3)

Y1 = rnorm(n/2, 5 + 0.25 * sin(X1 * 3 * pi), sd = sigma)
Y2 = rnorm(n, 3 + 0.5 * sin(X2 * 5 * pi), sd = sigma)
Y3 = rnorm(n/2, 1 + 0.25 * sin(X3 * 3 * pi), sd = sigma)
Y = c(Y1, Y2, Y3)

#scaled Y
Y.s = Y/(sd(Y)/sd(X))

#smoothing parameter (non-optimized)
h0 = 0.075

# Xgrid = seq(from = min(X), to = max(X), length.out = 100)
# Ygrid = seq(from = min(Y.s), to = max(Y.s), length.out = 5)
# mesh = expand.grid(Xgrid, Ygrid)
mesh = as.matrix(expand.grid(seq(from=0,to=1, length.out=100), c(1,2,3)))
#modal regression
RM_Y = modalReg(Xdata = X, Ydata = Y.s, hx = h0, hy = h0)
uni_Y = unimode(X, Y.s, hx = 0.1, hy = 0.1)
alpha=.95
RM_ps = quantile(abs(RM_Y-Y.s), alpha)
uni_ps = quantile(abs(uni_Y-Y.s), alpha)

pdf("./Figures/unimodalVSmultimodal.pdf", height = 7,
    width = 12)
par(mfrow=c(1,2))
plot(X,Y.s, cex = 0.4, ylab="Y", ylim=c(-0.1,1.2))
ordering=order(X)
polygon(x = c(X[ordering],rev(X[ordering])), 
        y = c((uni_Y-uni_ps)[ordering],rev((uni_Y+uni_ps)[ordering])),
        col = "orange", density = 50)
points(x = X, y = uni_Y, col = "dodgerblue", pch = 19, cex = 0.4)
points(x = X, y = uni_Y+uni_ps, col = "magenta", pch = 19, cex = 0.3)
points(x = X, y = uni_Y-uni_ps, col = "magenta", pch = 19, cex = 0.3)
legend("bottomleft", c("Unimodal Regression",
                       paste(100*alpha, "% PS", sep="")),
       lwd = c(2,2), col=c("dodgerblue","magenta"), cex = 0.6, lty = c(1,2) )


plot(X,Y.s, cex = 0.4, ylab="Y", ylim=c(-0.1,1.2))
order1=order(X[1:250])
order2=order(X[251:750])
order3=order(X[751:1000])
polygon(x = c(X[order1],rev(X[order1])), 
        y = c((RM_Y-RM_ps)[order1],
              rev((RM_Y+RM_ps)[order1])),
        col = "pink", density = 50)
polygon(x = c(X[251:750][order2], rev(X[251:750][order2])), 
        y = c((RM_Y-RM_ps)[251:750][order2],
              rev((RM_Y+RM_ps)[251:750][order2])),
        col = "pink", density = 50)
polygon(x = c(X[751:1000][order3], rev(X[751:1000][order3])), 
        y = c((RM_Y-RM_ps)[751:1000][order3],
              rev((RM_Y+RM_ps)[751:1000][order3])),
        col = "pink", density = 50)
points(x=X, y=RM_Y, col="blue", pch=19, cex = 0.4)
points(x = X, y = RM_Y+RM_ps, col = "red", pch = 19, cex = 0.3)
points(x = X, y = RM_Y-RM_ps, col = "red", pch = 19, cex = 0.3)
legend("bottomleft", c("Multimodal Regression",
                       paste(100*alpha, "% PS", sep="")),
       lwd = c(2,2), col=c("blue","red"), cex = 0.6, lty = c(1,2) )

dev.off()

