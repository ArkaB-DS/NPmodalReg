source("./R/meanshift1d.R")

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

Xgrid = seq(from = min(X), to = max(X), length.out = 100)
Ygrid = seq(from = min(Y.s), to = max(Y.s), length.out = 5)
mesh = expand.grid(Xgrid, Ygrid)

#modal regression
y.modalreg = modalReg(Xdata = X, Ydata = Y.s, xgrid = mesh[, 1],
                      ygrid = mesh[, 2], hx = h0, hy = h0)

pdf("./Figures/modmani.pdf", height = , width = 7)
plot(X, Y.s, cex =0.5,  ylim = c(-0.5, 1), ylab="Y")
lines(mesh[, 1][1:100], y.modalreg[1:100], lwd = 2, col = 'blue')
lines(mesh[, 1][401:500], y.modalreg[401:500], lwd = 2, col = 'blue',)
#locator(2)
text(0.05701435,0.7051011,"S1", col="red",cex=1.5)
text(0.05701435,0.2200192,"S2", col="red",cex=1.5)
dev.off()

