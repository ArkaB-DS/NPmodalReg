unimode<-function(Xdata, Ydata, xgrid=Xdata, ygrid=Ydata, hx=0.01, hy=0.01)
{
  n <- length(Xdata)
  mat <- matrix(0, nrow = length(xgrid), ncol = length(ygrid))
  for(i in 1:length(xgrid))
  {
    for(j in 1:length(ygrid))
    {
      mat[i,j] = 1/(n*hx*hy)*sum(dnorm(Xdata-xgrid[i], sd = hx) * 
                                   dnorm(Ydata-ygrid[j], sd = hy))
    }
  }
  # mat <- outer(xgrid, ygrid,
  #              function(x, y) 1/(n*hx*hy)*sum(dnorm(Xdata-x, sd = hx) * 
  #                                               dnorm(Ydata-y, sd = hy)) )
  ygrid[apply(mat, 1, which.max)]
}
