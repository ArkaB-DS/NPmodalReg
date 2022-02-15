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

par(mfrow=c(1,2))
plot(X,Y.s, cex = 0.4, ylab="Y")
points(x = X, y = uni_Y, col = "red", pch = 19, cex = 0.6)
plot(X,Y.s, cex = 0.4, ylab="Y")
points(x=X, y=RM_Y, col="blue", pch=19, cex = 0.5)


