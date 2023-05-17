# R version: 4.0.3 (2020-10-10)
# Platform: x86_64-apple-darwin17.0
# OS: macOS Mojave 10.14.1

#########################################################################
### Brownian Motion Simulation

set.seed(72)


# t is the endpoint of our time range,
# n is the number of steps in the simulation
# So in practice we're sampling n copies of the normal distribution
# on a timescale of length t/n and summing them

t = 100
n = 1000
numPaths = 4
steps = seq(0, t, length = n + 1)
A = replicate(numPaths, {
  bm <- c(0, cumsum(rnorm(n, 0, sqrt(t/n))))
}) 

timescale = seq(0, t, length = n + 1)

plot(timescale, A[, 1], type = "l", col = "blue", ylim = c(-20, 20), xlab = "t", ylab = "W(t)")
lines(timescale, A[, 2], type = "l", col = "orange")
lines(timescale, A[, 3], type = "l", col = "red")
lines(timescale, A[, 4], type = "l", col = "green")
