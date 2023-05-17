# R version: 4.0.3 (2020-10-10)
# Platform: x86_64-apple-darwin17.0
# OS: macOS Mojave 10.14.1



#########################################################################
### Toy simulation for testing the Nelson-Aalen estimator and the one
# and two sample log-rank tests



# Set rng seed for predictable results
set.seed(72)



#########################################################################
### Constants & Functions



# Gompertz coefficients for the two populations
b1 = 0.00001088387
c1 = 0.1043143
b2 = 0.000001619257
c2 = 0.1213714




# Forces of mortality
fom1 = function(t) {
  return(b1 * exp(c1 * t))
}

fom2 = function(t) {
  return(b2 * exp(c2 * t))
}



# Cumulative intensities
Lambda1 = function(t) {
  return((b1/c1) * (exp(c1 * t) - 1))
}

Lambda2 = function(t) {
  return((b2/c2) * (exp(c2 * t) - 1))
}



# Inverse distribution functions
Inverse1 = function(y) {
  return(log(1 - ((c1/b1) * log(1 - y))) / c1)
}

Inverse2 = function(y) {
  return(log(1 - ((c2/b2) * log(1 - y))) / c2)
}



#########################################################################
### Estimators & Tests



### Nelson-Aalen estimator
# We assume the sample is a vector of lifetimes,
# sorted in ascending order
NelsonAalen = function(t, smpl) {
  A = 0
  Y = length(smpl)
  
  # We use the following scheme throughout the code to
  # loop over jump times of the N corresponding to smpl up to time t
  # The if will hit exactly as many times
  # as there are ended lifetimes i.e. jumps before t
  # Because smpl is sorted, the exact time for the i:th jump
  # is given by smpl[i] and we can break once a lifetime ecxeeds t
  # This is not particularly efficient (you really should at least
  # be caching some stuff) but it gets the job done
  for (i in 1:Y) {
    if (smpl[i] <= t) {
      A = A + (1 / (Y - i + 1))
    } else {
      break()
    }
  }
  return(A)
}



### One sample log-rank test
# The statistic Z(t) is N(t) - the integral term, let us call it E(t)
# We calculate these separately

# Simply the number of deaths by time t
Nt = function(t, smpl) {
  N = 0
  for (i in 1:length(smpl)) {
    if (smpl[i] <= t) {
      N = i
    } else {
      break()
    }
  }
  return(N)
}

# Number at risk at time t
atRisk = function(t, smpl) {
  Y = length(smpl)
  for (i in 1:length(smpl)) {
    if (smpl[i] <= t) {
      Y = Y - 1
    } else {
      break()
    }
  }
  return(Y)
}

# Integral of Y(s)*alpha(s)
E = function(t, smpl, intensity) {
  integrand = function(s) {
    return(atRisk(s, smpl) * intensity(s))
  }
  return(integrate(Vectorize(integrand), lower = 0, upper = t, subdivisions = 2000, rel.tol=.Machine$double.eps^.05)$value)
}

# For the log-rank choice of K, the variance estimate reduces to
# E(t) so we may now calculate the test statistic:
logRank1 = function(t, smpl, intensity) {
  Et = E(t, smpl, intensity)
  return((Nt(t, smpl) - Et) / sqrt(Et))
}



### Two sample log-rank test
# The Z(t) is now a difference of two
# counting process integral functions featuring Y1 and Y2
# For convinience we assume the samples to be of the same size
Z2 = function(t, smpl1, smpl2) {
  # First term, an integral by N1
  s1 = 0
  for (i in 1:length(smpl1)) {
    if (smpl1[i] < t) {
      Y1 = atRisk(smpl1[i], smpl1)
      Y2 = atRisk(smpl1[i], smpl2)
      if (Y1 + Y2 > 0) {
        s1 = s1 + (Y2 / (Y1 + Y2))
      }
    } else {
      break()
    }
  }
  
  # Second term, an integral by N2
  s2 = 0
  for (i in 1:length(smpl2)) {
    if (smpl2[i] < t) {
      Y1 = atRisk(smpl2[i], smpl1)
      Y2 = atRisk(smpl2[i], smpl2)
      if (Y1 + Y2 > 0) {
        s2 = s2 + (Y1 / (Y1 + Y2))
      }
    } else {
      break()
    }
  }
  
  return(s1 - s2)
}

# Variance estimate
sigma2 = function(t, smpl1, smpl2) {
  # This integral is wrt the sum of N1 and N2
  # We rearrange the terms in the sum so that we first loop over
  # all jump times of N1 and then those of N2
  s = 0
  for (i in 1:length(smpl1)) {
    if (smpl1[i] < t) {
      Y1 = atRisk(smpl1[i], smpl1)
      Y2 = atRisk(smpl1[i], smpl2)
      if (Y1 + Y2 > 0) {
        s = s + ((Y1 * Y2) / ((Y1 + Y2)**2))
      }
    } else {
      break()
    }
  }
  
  for (i in 1:length(smpl2)) {
    if (smpl2[i] < t) {
      Y1 = atRisk(smpl2[i], smpl1)
      Y2 = atRisk(smpl2[i], smpl2)
      if (Y1 + Y2 > 0) {
        s = s + ((Y1 * Y2) / ((Y1 + Y2)**2))
      }
    } else {
      break()
    }
  }
  
  return(s)
}

# The final statistic
logRank2 = function(t, smpl1, smpl2) {
  return(Z2(t, smpl1, smpl2) / sqrt(sigma2(t, smpl1, smpl2)))
}

#########################################################################
### Simulation



# Sample size
n = 5000



# Vectors of lifetimes
# Samples 1 and 2 will be from population 1, sample 3 from population 2
sample1 = rep(0, n)
sample2 = rep(0, n)
sample3 = rep(0, n)



# Simulate the lifetimes
# (R would also let you do this more "elegantly" without the loop
# but I find this more intuitively readable)
for (i in 1:n) {
  sample1[i] = Inverse1(runif(1))
  sample2[i] = Inverse1(runif(1))
  sample3[i] = Inverse2(runif(1))
}



# Histograms
hist(sample1, breaks = 30, xlim = c(0, 120), ylim = c(0, 1200))
hist(sample2, breaks = 30, xlim = c(0, 120), ylim = c(0, 1200))
hist(sample3, breaks = 30, xlim = c(0, 120), ylim = c(0, 1200))



# Sort the samples for easier computation
sample1 = sort(sample1)
sample2 = sort(sample2)
sample3 = sort(sample3)



# Calculate Nelson-Aalen estimates
# We use a simple "yearly" evaluation density,
# i.e. we calculate and store the values of the estimator
# at t = 1, 2, ..., 115
endAge = 115
timescale = 1:endAge
A1 = rep(0, length(timescale))
A2 = rep(0, length(timescale))
A3 = rep(0, length(timescale))

# (The i is redundant with this configuration but
# would be needed for different timescales)
i = 1
for (t in timescale) {
  A1[i] = NelsonAalen(t, sample1)
  
  
  A2[i] = NelsonAalen(t, sample2)
  
  
  A3[i] = NelsonAalen(t, sample3)
  
  
  i = i + 1
}



#########################################################################
### Results



# Plot these against the actual cumulative intensities
title = "Nelson-Aalen-estimaatit kumulatiiviselle intensiteetille otoksissa 1 ja 3"
xrange = timescale
plot(timescale, A1, type = "s", xlab = "Age", ylab = "A(t)", main = title, col = "blue")
lines(timescale, A3, type = "s", col = "orange")
lines(Lambda1(xrange), type = "l", col = "black")
lines(Lambda2(xrange), type = "l", col = "gray")
legend(x = "topleft", legend=c("Otos 1", "Otos 3", "Lambda 1", "Lambda 2"), fill = c("blue", "orange", "black", "gray"))



# One sample log-rank test for sample 1 versus its actual intensity
lr1 = logRank1(endAge, sample1, fom1)
print(lr1)
print(2 * (1 - pnorm(abs(lr1))))

# p = 0.8211819, so we correctly fail to reject the null hypothesis

# One sample log-rank test for sample 1 versus
# the false intensity of population 2
lr1 = logRank1(endAge, sample1, fom2)
print(lr1)
print(2 * (1 - pnorm(abs(lr1))))

# p is so small it rounds out to zero,
# so the null is rejected as it should be

# Two sample log-rank test for sample 1 versus sample 2
lr2 = logRank2(endAge, sample1, sample2)
print(lr2)
print(2 * (1 - pnorm(abs(lr2))))

# p = 0.3463391, null correctly holds

# Two sample log-rank test for sample 1 versus sample 3
lr2 = logRank2(endAge, sample1, sample3)
print(lr2)
print(2 * (1 - pnorm(abs(lr2))))

# again p is approximately 0 and null is correctly rejected


