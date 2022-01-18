## Section 2.2 Histogram
days <- c(1, 27, 49, 84, 126, 257, 1, 27, 49, 84, 129, 311,
	5, 30, 54, 84, 134, 314, 7, 30, 56, 90, 144, 322,
	8, 31, 56, 91, 147, 369, 8, 31, 62, 92, 153, 415,
	13, 32, 63, 93, 163, 573, 14, 34, 65, 93, 167, 609,
	14, 35, 65, 103, 175, 640, 17, 36, 67, 103, 228, 737,
	18, 37, 75, 111, 231, 21, 38, 76, 112, 235,
	21, 39, 79, 119, 242, 22, 39, 82, 122, 256)

# Fig 2.2
hist(days, freq = FALSE, xlim = c(-200, 1000), ylim = c(0, 0.01),
main = "Histogram of lengths of treatment of control patients in suicide study.",
xlab = "Length of treatment (days)", ylab = "Relatve frequency",
breaks = 13)

# Fig 2.1
hist(faithful$eruptions, freq = FALSE, xlim = c(0,6), ylim = c(0,0.6),
main = "Histogram of eruption lengths of Old Faithful geyser.",
xlab = "Eruption length (min)", ylab = "Relatve frequency")

## Section 2.3 Naive Estimator

naiveDE <- function (x, data, bin.width = 0.25)
{
  return (sapply(x, function (x) 1/length(data)*
                  sum( 1/bin.width*0.5*(abs( (x-data) / bin.width ) <=1) ) ) )
}

x = seq(min(faithful$eruptions) - 0.5, max(faithful$eruptions) + 0.5, by = 1e-5)
y = naiveDE(x , faithful$eruptions)

# Fig 2.3
plot(x, y, type = "l", xlim= c(0,6), ylim = c(0, 0.7), xlab = "Eruption length (min)",
  ylab = "Density estimate", main = "Naive estimate constructed from Old Faithful geyser data,
  h=0.25")

# Section 2.4 Kernel Estimator

kernelDE <- function(x, data, bandwidth = 0.25, kernel)
{
  return (sapply(x, function (x) 1/length(data)*
                   sum( 1/bandwidth*kernel((x-data)/bandwidth) ) ) )
}

x = seq(min(faithful$eruptions) - 0.5, max(faithful$eruptions) + 0.5, by = 1e-5)
y = kernelDE(x , faithful$eruptions, kernel = function (k) dnorm(k))

# Fig 2.8
plot(x, y, type = "l", xlim= c(0,6), ylim = c(0, 0.6), xlab = "Eruption length (min)",
     ylab = "Density estimate", main = "Kernel estimate for Old Faithful geyser data,
  window width 0.25")

x = seq(-200, 1000, by = 1e-3)

# Fig 2.9
par(mfcol = c(2,1))

y = kernelDE(x, days, 20, kernel = function (k) dnorm(k))
plot(x, y, type = "l", ylim= c(0,0.007), xlim = c(-200, 1000), xlab = "Length of treatment (days)",
     ylab = "Density estimate")
title(main = "Kernel estimate for suicide study data. \nWindow widths: (a) 20; (b) 60.")
y = kernelDE(x, days, 60, kernel = function (k) dnorm(k))
plot(x, y, type = "l", ylim= c(0,0.007), xlim = c(-200, 1000), xlab = "Length of treatment (days)",
     ylab = "Density estimate")



