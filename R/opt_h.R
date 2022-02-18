source("./R/meanshift1d.R")

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


### Section 3: Bandwidth selection via minimizing prediction set
alpha=0.95

x_n = 50	
y_n = 20

x_gr =seq(from=0,to=1, length.out=x_n)
y_gr = seq(from=0.5, to=3.5, length.out=y_n)
# mesh points

k_RS = rep(NA, x_n)
# number of clusters

h_seq = seq(from=0.02, to=0.20, length.out = 100)
# test a series of h: from 0.02 to 0.20

RM_list = list()
PI_seq = rep(NA, length(h_seq))

for(i_tmp in 1:length(h_seq)){
  h_x =h_seq[i_tmp]
  h_y = h_x*sd(Y)/sd(X)
  
  RM_Y = modalReg(X,Y,X,Y,h_x,h_y)
  RM_list[[i_tmp]] = RM_Y
  
  
  for(i in 1:x_n){
    pt_tmp = cbind(rep(x_gr[i],y_n), y_gr)
    # initialize points at each x grid
    
    rm_y_tmp = modalReg(X,Y,rep(x_gr[i],y_n),y_gr, h_x, h_y)	
    # finding conditional local modes
    
    y_tmp = round(rm_y_tmp, digits=5)
    k_RS[i] = length(unique(y_tmp))
    # number of local mode at the i-th mesh point
  }
  
  q_RS = quantile(abs(RM_Y-Y), alpha)
  # width of prediction set
  
  PI_seq[i_tmp] = mean(k_RS)*q_RS
  # area of prediction set
}

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

pdf("./Figures/optH.pdf", height=7, width=7)
plot(x=h_seq, y= PI_seq, lwd=2, type="l",
     ylab=paste("Size of ",100*alpha,"% Prediction Sets", sep=""),
     xlab="h", col = "blue")
abline(v= h_seq[which(PI_seq==min(PI_seq))], col="red", lwd=2, lty=2)
abline(h= quantile(abs(loess_fit$res),alpha), col="black", lwd=2, lty=3)
legend("topright", c("Modal Regression","Local Regression", "Optimal h"),
       col=c("blue","black","red"), lwd=c(2,2,2), lty=c(1,3,2),
       cex=0.8)
dev.off()

# plot(x=h_seq, y= PI_seq, lwd = 2, type="l",
#      ylab="Size of prediction set", xlab="h",col='blue',
#      main=paste("Size of ",100*alpha,"% Prediction interval", sep=""))
# abline(v= h_seq[which(PI_seq==min(PI_seq))], col="red", lty=2)
# 
# 

