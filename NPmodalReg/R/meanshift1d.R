modalReg <- function(Xdata, Ydata, xgrid, ygrid, hx, hy, tol = 1e-6, iter = 100)
{
  # initialization
  x0 = xgrid
  y0 = ygrid
  
  i = 0
  diff = 0.5
  
  while(i < iter & diff > tol)
  {
    ytemp = y0
    y0 = sapply(1:length(y0), 
                function(x) sum(Ydata*dnorm(y0[x] - Ydata, sd = hy) * 
                                  dnorm(x0[x] - Xdata, sd = hx)) / 
                  sum(dnorm(y0[x] - Ydata, sd = hy) * dnorm(x0[x] - Xdata, sd = hx)))
  }
  i = i + 1
  diff = max(abs(ytemp - y0))
  
  return(y0)
}
