# Comparison of Modal Regression with Local Regression

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

MS_fit = modalReg(X, Y.s, hx = h0, hy = h0)
#predicted Y on the modal curves of each data point

alpha=0.95

MS_ps = quantile(abs(MS_fit-Y.s), alpha)
# the quantile of loss; this is larger than the Hausdorff distance so that the prediction set has proper coverage

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

plot(X, Y.s, cex = 0.5, ylim = c(-1, 1), ylab="Y")
lines(X[order(X)], fit_loc$fitted[order(X)],
      lwd=2, col="red")
lines(X[order(X)], fit_loc$fitted[order(X)]+quantile(fit_loc$res,alpha),
      lwd=2, col="blue",lty=2)
lines(X[order(X)], fit_loc$fitted[order(X)]-quantile(fit_loc$res,alpha),
      lwd=2, col="blue",lty=2)
legend("bottomleft", c("Local Regression",paste(100*alpha,"% PS, Local", sep="")),
       lwd=c(2,2), col=c("red","blue") , cex = 0.6, lty = c(1,2))

plot(X, Y.s, cex =0.5,  ylim = c(-1, 1), ylab="Y")
lines(mesh[, 1][1:100], y.modalreg[1:100], lwd = 2, col = 'red')
lines(mesh[, 1][401:500], y.modalreg[401:500], lwd = 2, col = 'red',)
lines(mesh[, 1][1:100], y.modalreg[1:100] + MS_ps, lwd = 2, col = "blue",lty=2)
lines(mesh[, 1][401:500], y.modalreg[401:500] + MS_ps, lwd=2, col = "blue",lty=2)
lines(mesh[, 1][1:100], y.modalreg[1:100] - MS_ps, lwd = 2, col = "blue",lty=2)
lines(mesh[, 1][401:500], y.modalreg[401:500] - MS_ps, lwd = 2, col = "blue",lty=2)
legend("bottomleft", c("Modal Regression",
                       paste(100*alpha, "% PS, Modal", sep="")),
       lwd = c(2,2), col=c("blue","dodgerblue"), cex = 0.6, lty = c(1,2) )


### Section 1: Generate Data
n1=100
n2=100
n3=100

sd_y = 0.2

x1 = runif(n1,0,0.4)
x2 = runif(n2,0.3, 0.7)
x3 = runif(n3,0.6,1.0)


y1 = rnorm(n1,3,sd=sd_y)
y2 = rnorm(n2,2,sd=sd_y)
y3 = rnorm(n3,1,sd=sd_y)


X = c(x1,x2,x3)
Y = c(y1,y2,y3)


### Section 2: Modal Regression
h_x = 0.1090909	#optimal choice
h_y = h_x*sd(Y)/sd(X)

Grids = as.matrix(expand.grid(seq(from=0,to=1, length.out=100), c(1,2,3)))

RM_Y = modalReg(X,Y,Grids[,1], Grids[,2],   h_x,h_y)


###
h_tmp = hclust(dist(cbind(Grids[,1],RM_Y)))
# hierachical clustering over the modal regression points
lab_tmp = cutree(h_tmp,h=1)
# clustering

clusters = list()
for( i in 1:max(lab_tmp)){
  w_tmp = which(lab_tmp==i)
  clusters[[i]] = cbind(Grids[w_tmp,1], RM_Y[w_tmp])
}
# each list element is the modal points for that cluster

alpha=0.95
# q_RS_opt = quantile(abs(RM_list[[which(PI_seq==min(PI_seq))]]-Y), alpha)
q_RS_opt = 0.4154116 


#### Section 6: Local regression

span_seq = seq(from=0.1, to=0.9, length.out=100)
loc_PS_seq = rep(NA,100)

for(i in 1:100){
  fit_loc_tmp = loess(Y~X,span=span_seq[i])
  loc_PS_seq[i] = quantile(fit_loc_tmp$res, alpha)
}
s_opt = span_seq[which(loc_PS_seq==min(loc_PS_seq))]
#the optimal one

loess_fit  = loess(Y~X,span=s_opt)

plot(X,Y,cex=0.5, ylim=c(0,4))
lines(loess_fit$x[order(loess_fit$x)],loess_fit$fitted[order(loess_fit$x)], 
      col="red", lwd=2)
lines(loess_fit$x[order(loess_fit$x)],loess_fit$fitted[order(loess_fit$x)]+quantile(abs(loess_fit$res),alpha),
      col="blue", lwd=2,lty=2)
lines(loess_fit$x[order(loess_fit$x)],loess_fit$fitted[order(loess_fit$x)]-quantile(abs(loess_fit$res),alpha),
      col="blue", lwd=2, lty=2)
legend("bottomleft",c("Local Regression", paste(alpha,"% PI", sep="")),
       lwd=c(2,2), col=c("red","blue"), cex=0.6, lty=c(1,2) )

plot(X,Y, cex=0.5,ylim=c(0,4))
for(i in 1:max(lab_tmp)){
  lines(clusters[[i]][order(clusters[[i]][,1]),], col="red",
        lwd=2)
  lines(x=clusters[[i]][order(clusters[[i]][,1]),][,1], y=clusters[[i]][order(clusters[[i]][,1]),][,2]+q_RS_opt,
        col="blue", lwd=2,lty=2)
  lines(x=clusters[[i]][order(clusters[[i]][,1]),][,1], y=clusters[[i]][order(clusters[[i]][,1]),][,2]-q_RS_opt,
        col="blue", lwd=2, lty=2)
}
legend("bottomleft",c("Modal Regression", paste(alpha,"% PI", sep="")),
       lwd=c(2,2), col=c("red","blue"), cex=0.6, lty=c(1,2) )

