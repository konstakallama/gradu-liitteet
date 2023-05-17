# R version: 4.0.3 (2020-10-10)
# Platform: x86_64-apple-darwin17.0
# OS: macOS Mojave 10.14.1



#########################################################################
### Data

# This is obviously my local wd, use your own
setwd("~/Documents/School/HY/MT/d")

# The data format is as follows:
# Both files are ordinary .csv:s,
# where the first row and column are labels.
# Rows go 1898 to 2001, these represent birth cohorts.
# Columns go form 2007 to 2020, these are the observation years.
# In the atRisk files, the value represents
# the number of people of that cohort alive at the start of that year.
# In the deaths files, the value represents
# the number of people of that cohort who died during that year.

YELatRiskDataM = read.csv("YELatRiskM.csv")
YELatRiskDataF = read.csv("YELatRiskF.csv")

TyELatRiskDataM = read.csv("TyELatRiskM.csv")
TyELatRiskDataF = read.csv("TyELatRiskF.csv")

YELdeathsDataM = read.csv("YELdeathsM.csv")
YELdeathsDataF = read.csv("YELdeathsF.csv")

TyELdeathsDataM = read.csv("TyELdeathsM.csv")
TyELdeathsDataF = read.csv("TyELdeathsF.csv")

dataStartYear = 2007
dataEndYear = 2020

dataStartCohort = 1898
dataEndCohort = 2001



#########################################################################
### Functions

# Constants needed for calculating the theorethical
# Force of Mortality for different cohorts/sexes.
# See Lehtovirta (2020), section 2.2 for the model.

b2 = function(cohort) {
  if (cohort < 1930) {
    return(5)
  } else if (cohort < 1940) {
    return(3)
  } else if (cohort < 1950) {
    return(2)
  } else if (cohort < 1960) {
    return(0)
  } else if (cohort < 1970) {
    return(-2)
  } else if (cohort < 1980) {
    return(-3)
  } else if (cohort < 1990) {
    return(-5)
  } else if (cohort < 2000) {
    return(-7)
  } else if (cohort < 2010) {
    return(-8)
  } else {
    return(-10)
  }
}

a1m = function(age, cohort) {
  if (age + b2(cohort) <= 70) {
    return(exp(((6/7) * 1.027) - 11.18))
  } else {
    return(exp(((6/7) * 1.217) - 12.68))
  }
}

a1f = function(age, cohort) {
  if (age + b2(cohort) <= 70) {
    return(exp(((6/7) * 1.031) - 11.86))
  } else {
    return(exp(((6/7) * 1.416) - 14.79))
  }
}

a2m = function(age, cohort) {
  if (age + b2(cohort) <= 70) {
    return((6/7) * 0.1027)
  } else {
    return((6/7) * 0.1217)
  }
}

a2f = function(age, cohort) {
  if (age + b2(cohort) <= 70) {
    return((6/7) * 0.1031)
  } else {
    return((6/7) * 0.1416)
  }
}



#########################################################################
# The theorethical forces of mortality

FoMm = function(age, cohort) {
  return(a1m(age, cohort) * exp(a2m(age, cohort) * (age + b2(cohort))))
}

FoMf = function(age, cohort) {
  return(a1f(age, cohort) * exp(a2f(age, cohort) * (age + b2(cohort))))
}



#########################################################################
# Utility functions for reading the data

atRisk = function(age, cohort, atRiskData) {
  return(atRiskData[cohort - dataStartCohort + 1, cohort + age - dataStartYear + 2])
}

dead = function(age, cohort, deathsData) {
  return(deathsData[cohort - dataStartCohort + 1, cohort + age - dataStartYear + 2])
}



#########################################################################
# Function for calculating the Nelson-Aalen estimator
# This gives an array of yearly estimator values
# for all the ages of the given cohorts that are present in the data
# For earlier/later ages the estimator will be
# 0/constant at its max value since there is no data to estimate from

NelsonAalenTieCorrected = function(startCohort, endCohort, atRiskData, deathsData) {
  minAge = max(dataStartYear - endCohort, 1)
  maxAge = max(dataEndYear - startCohort, 1)
  A = rep(0, maxAge - minAge + 1)
  
  for (age in minAge:maxAge) {
    d = 0
    Y = 0
    for (cohort in startCohort:endCohort) {
      if (age + cohort > dataStartYear - 1 && age + cohort < dataEndYear + 1) {
        d = d + dead(age, cohort, deathsData)
        Y = Y + atRisk(age, cohort, atRiskData)
      }
    }
    
    # Tie correction
    increment = 0
    if (d > 0) {
      for (i in 0:(d - 1)) {
        if (Y - i > 0) {
          increment = increment + (1 / (Y - i))
        }
      }
    }
    
    if (age == minAge) {
      A[age - minAge + 1] = increment
    } else {
      A[age - minAge + 1] = A[age - minAge] + increment
    }
  }
  
  return(A)
}



#########################################################################
# Analytical integrals of the forces of mortality
# These assume a two-part Gompertz model
# as specified in Lehtovirta (2020)

LambdaM = function(t, cohort) {
  b = b2(cohort)
  a11 = a1m(0, cohort)
  a21 = a2m(0, cohort)
  a12 = a1m(100, cohort)
  a22 = a2m(100, cohort)
  
  if (t <= 70 - b) {
    return((a11/a21) * exp(a21 * b) * (exp(a21 * t) - 1))
  } else {
    return(((a11/a21) * exp(a21 * b) * (exp(a21 * (70 - b)) - 1)) + ((a12/a22) * exp(a22 * b) * (exp(a22 * t) - exp(a22 * (70 - b)))))
  }
}

LambdaF = function(t, cohort) {
  b = b2(cohort)
  a11 = a1f(0, cohort)
  a21 = a2f(0, cohort)
  a12 = a1f(100, cohort)
  a22 = a2f(100, cohort)
  
  if (t <= 70 - b) {
    return((a11/a21) * exp(a21 * b) * (exp(a21 * t) - 1))
  } else {
    return(((a11/a21) * exp(a21 * b) * (exp(a21 * (70 - b)) - 1)) + ((a12/a22) * exp(a22 * b) * (exp(a22 * t) - exp(a22 * (70 - b)))))
  }
}



#########################################################################
# Functions for calculating the one-sample log-rank statistic
# Given that the force of mortality in the model
# changes with every birth decade,
# the easiest way this makes sense is to assume that the
# cohort range given here as a parameter falls inside the same decade.
# In practice this just uses the force of mortality for
# the first cohort so it won't break if you give it a longer range,
# but it will overestimate the force of mortality
# and may give results that are hard to interpret.

# By default these functions will give
# the value of the statistic at the largest possible time,
# that is the age of the first cohort in the last observation year.

# First the integral of Y(s)*alpha(s) aka E:

E = function(startCohort, endCohort, atRiskData, deathsData, integrand) {
  startAge = max(dataStartYear - endCohort, 1)
  endAge = max(dataEndYear - startCohort, 1)
  t = endAge
  
  s = 0
  for (age in startAge:t) {
    Y = 0
    for (cohort in startCohort:endCohort) {
      if (cohort + age >= dataStartYear && cohort + age <= dataEndYear) {
        # The atRisk numbers in the data correspond to
        # the population size at the start of the year
        # We assume the deaths are spread evenly throughout the year, 
        # so that the value of Y at time s is approximately
        # the midpoint between the risk set sizes at the start
        # and end of the year
        # We use the deaths data for this, because we also
        # assume that censorings happen only at the start of the year
        # (as they do in our case where censorings are
        # due to study time running out)
        Y = Y + (((2 * atRisk(age, cohort, atRiskData)) - dead(age, cohort, deathsData)) / 2)
      }
    }
    intv = integrate(Vectorize(integrand), lower = age, upper = age + 1)$value
    s = s + (Y * intv)
  }
  return(s)
}

Em = function(startCohort, endCohort, atRiskData, deathsData, fomCohort) {
  integrand = function(x) {
    return(FoMm(x, fomCohort))
  }
  return(E(startCohort, endCohort, atRiskData, deathsData, integrand))
}

Ef = function(startCohort, endCohort, atRiskData, deathsData, fomCohort) {
  integrand = function(x) {
    return(FoMf(x, fomCohort))
  }
  return(E(startCohort, endCohort, atRiskData, deathsData, integrand))
}

# N aka simply the number of deaths by age:

N = function(startCohort, endCohort, deathsData) {
  startAge = max(dataStartYear - endCohort, 1)
  endAge = max(dataEndYear - startCohort, 1)
  t = endAge
  
  s = 0
  for (age in startAge:t) {
    for (cohort in startCohort:endCohort) {
      if (cohort + age >= dataStartYear && cohort + age <= dataEndYear) {
        s = s + dead(age, cohort, deathsData)
      }
    }
  }
  return(s)
}

LogRankM = function(startCohort, endCohort, atRiskData, deathsData, fomCohort) {
  Nt = N(startCohort, endCohort, deathsData)
  Et = Em(startCohort, endCohort, atRiskData, deathsData, fomCohort)
  return((Nt - Et) / sqrt(Et))
}

LogRankF = function(startCohort, endCohort, atRiskData, deathsData, fomCohort) {
  Nt = N(startCohort, endCohort, deathsData)
  Et = Ef(startCohort, endCohort, atRiskData, deathsData, fomCohort)
  return((Nt - Et) / sqrt(Et))
}



#########################################################################

# Functions for the two-sample log-rank statistic
# Here we can use cohort ranges as long as we want
# since there's no changing force of mortality model to worry about

# Z = sum(dN1 * K/Y1) - sum(dN2 * K/Y2)
# K = Y1 * Y2 / Y.
# -> Z = sum(dN1 * Y2/Y.) - sum(dN2 * Y1/Y.)

Z2 = function(startCohort, endCohort, atRiskData1, atRiskData2, deathsData1, deathsData2) {
  startAge = max(dataStartYear - endCohort, 1)
  endAge = max(dataEndYear - startCohort, 1)
  t = endAge
  
  s = 0
  for (age in startAge:t) {
    for (cohort in startCohort:endCohort) {
      if (cohort + age >= dataStartYear && cohort + age <= dataEndYear) {
        N1 = dead(age, cohort, deathsData1)
        N2 = dead(age, cohort, deathsData2)
        Y1 = atRisk(age, cohort, atRiskData1)
        Y2 = atRisk(age, cohort, atRiskData2)
        
        # Similarly to the one-sample case, we approximate the size of
        # the "other" risk set at the jump times throughout the year
        # at the midpoint of its sizes at the start and end of the year
        Y1e = ((2*atRisk(age, cohort, atRiskData1)) - dead(age, cohort, deathsData1)) / 2
        Y2e = ((2*atRisk(age, cohort, atRiskData2)) - dead(age, cohort, deathsData2)) / 2
        
        s1 = 0
        if (Y1 + Y2e > 0) {
          s1 = ((N1 * Y2e) / (Y1 + Y2e))
        }
        
        s2 = 0
        if (Y1e + Y2 > 0) {
          s2 = ((N2 * Y1e) / (Y1e + Y2))
        }
        
        s = s + s1 - s2
      }
    }
  }
  return(s)
}

# variance estimate = sig^2 = sum(dN. * Y1 * Y2 / Y.^2)

sig2 = function(startCohort, endCohort, atRiskData1, atRiskData2, deathsData1, deathsData2) {
  startAge = max(dataStartYear - endCohort, 1)
  endAge = max(dataEndYear - startCohort, 1)
  t = endAge
  
  # Here we use doubles to prevent overflows
  s = 0.0
  for (age in startAge:t) {
    for (cohort in startCohort:endCohort) {
      if (cohort + age >= dataStartYear && cohort + age <= dataEndYear) {
        N1 = as.double(dead(age, cohort, deathsData1))
        N2 = as.double(dead(age, cohort, deathsData2))
        Y1 = as.double(atRisk(age, cohort, atRiskData1))
        Y2 = as.double(atRisk(age, cohort, atRiskData1))
        if (!(Y1 + Y2 == 0)) {
          s = s + (((N1 + N2) * Y1 * Y2) / ((Y1 + Y2)**2))
        }
      }
    }
  }
  return(s)
}

LogRank2 = function(startCohort, endCohort, atRiskData1, atRiskData2, deathsData1, deathsData2) {
  Zt = Z2(startCohort, endCohort, atRiskData1, atRiskData2, deathsData1, deathsData2)
  sigt = sig2(startCohort, endCohort, atRiskData1, atRiskData2, deathsData1, deathsData2)
  return(Zt / sqrt(sigt))
}



#########################################################################
#-----------------------------------------------------------------------#
#########################################################################
### Nelson-Aalen estimators


fullAgeRange = max(dataStartYear - dataEndCohort, 1):max(dataEndYear - dataStartCohort, 1)
naYELm = NelsonAalenTieCorrected(dataStartCohort, dataEndCohort, YELatRiskDataM, YELdeathsDataM)
naYELf = NelsonAalenTieCorrected(dataStartCohort, dataEndCohort, YELatRiskDataF, YELdeathsDataF)
naTyELm = NelsonAalenTieCorrected(dataStartCohort, dataEndCohort, TyELatRiskDataM, TyELdeathsDataM)
naTyELf = NelsonAalenTieCorrected(dataStartCohort, dataEndCohort, TyELatRiskDataF, TyELdeathsDataF)

# This code will produce plots just for viewing
# To save each plot to a .png, the following code was used:
# png(filename = paste(*filename*, ".png", sep = ""), width = 1200,
# height = 1200, units = "px", pointsize = 36, bg = "white",
# type = c("cairo-png"), symbolfamily="default")
# *plot code*
# dev.off()

par(mfrow = c(1, 1))

# Nelson-Aalen estimators for YEL vs TyEL males
plot(fullAgeRange, naYELm, type = "s", xlab = "Age", ylab = "A(t)", main = "", col = "blue", log = "")
lines(fullAgeRange, naTyELm, type = "s", col = "orange")
legend(x = "topleft", legend=c("TyEL", "YEL"), fill = c("orange", "blue"))

# Nelson-Aalen estimators for YEL vs TyEL males, logarithmic scale
plot(fullAgeRange, naYELm, type = "s", xlab = "Age", ylab = "A(t)", main = "", col = "blue", log = "y", ylim = c(0.0001, 10))
lines(fullAgeRange, naTyELm, type = "s", col = "orange")
legend(x = "topleft", legend=c("TyEL", "YEL"), fill = c("orange", "blue"))



# Nelson-Aalen estimators for YEL vs TyEL females
plot(fullAgeRange, naYELf, type = "s", xlab = "Age", ylab = "A(t)", main = "", col = "blue", log = "")
lines(fullAgeRange, naTyELf, type = "s", col = "orange")
legend(x = "topleft", legend=c("TyEL", "YEL"), fill = c("orange", "blue"))

# Nelson-Aalen estimators for YEL vs TyEL females, logarithmic scale
plot(fullAgeRange, naYELf, type = "s", xlab = "Age", ylab = "A(t)", main = "", col = "blue", log = "y", ylim = c(0.0001, 10))
lines(fullAgeRange, naTyELf, type = "s", col = "orange")
legend(x = "topleft", legend=c("TyEL", "YEL"), fill = c("orange", "blue"))



cohortPlotM = function(startCohort) {
  endCohort = startCohort + 9
  ageRange = (dataStartYear - endCohort):(dataEndYear - startCohort)
  
  naYELm = NelsonAalenTieCorrected(startCohort, endCohort, YELatRiskDataM, YELdeathsDataM)
  naTyELm = NelsonAalenTieCorrected(startCohort, endCohort, TyELatRiskDataM, TyELdeathsDataM)
  
  LmbdM = function(t) {
    return(LambdaM(t, startCohort) - LambdaM(ageRange[1] - 1, startCohort))
  }

  title = paste(startCohort, "-", endCohort, ", Miehet", sep = "")
  plot(ageRange, naYELm, type = "s", xlab = "Age", ylab = "A(t)", main = title, col = "blue", log = "")
  lines(ageRange, naTyELm, type = "s", col = "orange")
  lines(ageRange, Vectorize(LmbdM)(ageRange), type = "l", col = "black")
  legend(x = "topleft", legend=c("TyEL", "YEL", "Lambda"), fill = c("orange", "blue", "black"))
}

cohortPlotF = function(startCohort) {
  endCohort = startCohort + 9
  ageRange = (dataStartYear - endCohort):(dataEndYear - startCohort)
  
  naYELf = NelsonAalenTieCorrected(startCohort, endCohort, YELatRiskDataF, YELdeathsDataF)
  naTyELf = NelsonAalenTieCorrected(startCohort, endCohort, TyELatRiskDataF, TyELdeathsDataF)
  
  LmbdF = function(t) {
    return(LambdaF(t, startCohort) - LambdaF(ageRange[1] - 1, startCohort))
  }
  
  title = paste(startCohort, "-", endCohort, ", Naiset", sep = "")
  
  plot(ageRange, naYELf, type = "s", xlab = "Age", ylab = "A(t)", main = title, col = "blue", log = "")
  lines(ageRange, naTyELf, type = "s", col = "orange")
  lines(ageRange, Vectorize(LmbdF)(ageRange), type = "l", col = "black")
  legend(x = "topleft", legend=c("TyEL", "YEL", "Lambda"), fill = c("orange", "blue", "black"))
}

# Estimators vs the theoretical cumulative force of mortality
# for selected cohorts

cohortPlotM(1910)
cohortPlotF(1910)

cohortPlotM(1920)
cohortPlotF(1920)

cohortPlotM(1930)
cohortPlotF(1930)

cohortPlotM(1940)
cohortPlotF(1940)

cohortPlotM(1950)
cohortPlotF(1950)

cohortPlotM(1960)
cohortPlotF(1960)



#########################################################################
### Log-rank tests

# One-sample tests
# In order for these to have a force of mortality that's defined
# as the same for all samples,
# we restrict ourselves to one cohort at a time.

YELlr1M = rep(0, 6)
YELlr1F = rep(0, 6)
TyELlr1M = rep(0, 6)
TyELlr1F = rep(0, 6)

i = 1
for (cohort in seq(1910, 1960, 10)) {
  YELlr1M[i] = LogRankM(cohort, cohort + 9, YELatRiskDataM, YELdeathsDataM, cohort)
  YELlr1F[i] = LogRankF(cohort, cohort + 9, YELatRiskDataF, YELdeathsDataF, cohort)
  TyELlr1M[i] = LogRankM(cohort, cohort + 9, TyELatRiskDataM, TyELdeathsDataM, cohort)
  TyELlr1F[i] = LogRankF(cohort, cohort + 9, TyELatRiskDataF, TyELdeathsDataF, cohort)
  i = i + 1
}

pYELm = 2 * (1 - pnorm(abs(YELlr1M)))
pYELf = 2 * (1 - pnorm(abs(YELlr1F)))
pTyELm = 2 * (1 - pnorm(abs(TyELlr1M)))
pTyELf = 2 * (1 - pnorm(abs(TyELlr1F)))

# Test statistics for cohorts 1910-1919, 1920-1929, etc
YELlr1M

YELlr1F

TyELlr1M

TyELlr1F

# p-values
pYELm

pYELf

pTyELm

pTyELf



# Two-sample tests
# Here we can do the whole cohort range, but we'll also do
# the 10-year cohorts in the interest of completeness

lr2M = rep(0, 7)
lr2F = rep(0, 7)

i = 1
for (cohort in seq(1910, 1960, 10)) {
  lr2M[i] = LogRank2(cohort, cohort + 9, YELatRiskDataM, TyELatRiskDataM, YELdeathsDataM, TyELdeathsDataM)
  lr2F[i] = LogRank2(cohort, cohort + 9, YELatRiskDataF, TyELatRiskDataF, YELdeathsDataF, TyELdeathsDataF)
  i = i + 1
}

lr2M[7] = LogRank2(dataStartCohort, dataEndCohort, YELatRiskDataM, TyELatRiskDataM, YELdeathsDataM, TyELdeathsDataM)
lr2F[7] = LogRank2(dataStartCohort, dataEndCohort, YELatRiskDataF, TyELatRiskDataF, YELdeathsDataF, TyELdeathsDataF)

plr2M = 2 * (1 - pnorm(abs(lr2M)))
plr2F = 2 * (1 - pnorm(abs(lr2F)))

lr2M

lr2F

plr2M

plr2F



#########################################################################
#-----------------------------------------------------------------------#
#########################################################################
### Ikasiirto-adjusted data

# We repeat the important calculations of the previous section using
# data that's been adjusted to remove the need to account for b2
# e.g. for the 1940 cohort, b2 = 2 so they are treated by the model
# as 2 years older than they actually are
# so to bake this into the data, we have adjusted all
# of their birth years down by 2

# Read the new data

adjYELatRiskDataM = read.csv("adjYELatRiskM.csv")
adjYELatRiskDataF = read.csv("adjYELatRiskF.csv")

adjTyELatRiskDataM = read.csv("adjTyELatRiskM.csv")
adjTyELatRiskDataF = read.csv("adjTyELatRiskF.csv")

adjYELdeathsDataM = read.csv("adjYELdeathsM.csv")
adjYELdeathsDataF = read.csv("adjYELdeathsF.csv")

adjTyELdeathsDataM = read.csv("adjTyELdeathsM.csv")
adjTyELdeathsDataF = read.csv("adjTyELdeathsF.csv")

# Here we have new cohort ranges

dataStartCohort = 1893
dataEndCohort = 2011

#########################################################################
### Nelson-Aalen estimators for the adjusted data

# Now the youngest age for which we have any actual data is 10,
# so we'll use that as the starting point for the x-axis

adjFullAgeRange = 10:max(dataEndYear - dataStartCohort, 1)
adjnaYELm = NelsonAalenTieCorrected(dataStartCohort, dataEndCohort, adjYELatRiskDataM, adjYELdeathsDataM)
adjnaYELf = NelsonAalenTieCorrected(dataStartCohort, dataEndCohort, adjYELatRiskDataF, adjYELdeathsDataF)
adjnaTyELm = NelsonAalenTieCorrected(dataStartCohort, dataEndCohort, adjTyELatRiskDataM, adjTyELdeathsDataM)
adjnaTyELf = NelsonAalenTieCorrected(dataStartCohort, dataEndCohort, adjTyELatRiskDataF, adjTyELdeathsDataF)

# The 1950 cohort has b2 = 0, so we'll use that here
# What we've effectively done is normalize everyone so that they behave
# as if they were born in the 50s (as far as the model is concerned)

sLambdaM = function(t) {
  return(LambdaM(t, 1950) - LambdaM(adjFullAgeRange[1] - 1, 1950))
}

sLambdaF = function(t) {
  return(LambdaF(t, 1950) - LambdaF(adjFullAgeRange[1] - 1, 1950))
}



# Plots
par(mfrow = c(1, 1))

# YEL vs TyEL males
plot(adjFullAgeRange, adjnaYELm[adjFullAgeRange], type = "s", xlab = "Age", ylab = "A(t)", main = "", col = "blue", log = "")
lines(adjFullAgeRange, adjnaTyELm[adjFullAgeRange], type = "s", col = "orange")
lines(adjFullAgeRange, Vectorize(sLambdaM)(adjFullAgeRange), type = "l", col = "black")
legend(x = "topleft", legend=c("TyEL", "YEL", "Lambda"), fill = c("orange", "blue", "black"))

# YEL vs TyEL males, logarithmic scale
plot(adjFullAgeRange, adjnaYELm[adjFullAgeRange], type = "s", xlab = "Age", ylab = "A(t)", main = "", col = "blue", log = "y", ylim = c(0.0001, 10))
lines(adjFullAgeRange, adjnaTyELm[adjFullAgeRange], type = "s", col = "orange")
lines(adjFullAgeRange, Vectorize(sLambdaM)(adjFullAgeRange), type = "l", col = "black")
legend(x = "topleft", legend=c("TyEL", "YEL", "Lambda"), fill = c("orange", "blue", "black"))



# YEL vs TyEL females
plot(adjFullAgeRange, adjnaYELf[adjFullAgeRange], type = "s", xlab = "Age", ylab = "A(t)", main = "", col = "blue", log = "")
lines(adjFullAgeRange, adjnaTyELf[adjFullAgeRange], type = "s", col = "orange")
lines(adjFullAgeRange, Vectorize(sLambdaF)(adjFullAgeRange), type = "l", col = "black")
legend(x = "topleft", legend=c("TyEL", "YEL", "Lambda"), fill = c("orange", "blue", "black"))

# YEL vs TyEL females, logarithmic scale
plot(adjFullAgeRange, adjnaYELf[adjFullAgeRange], type = "s", xlab = "Age", ylab = "A(t)", main = "", col = "blue", log = "y", ylim = c(0.0001, 10))
lines(adjFullAgeRange, adjnaTyELf[adjFullAgeRange], type = "s", col = "orange")
lines(adjFullAgeRange, Vectorize(sLambdaF)(adjFullAgeRange), type = "l", col = "black")
legend(x = "topleft", legend=c("TyEL", "YEL", "Lambda"), fill = c("orange", "blue", "black"))



#########################################################################
### Tests for the adjusted data

# One-sample tests

adjYELlr1M = LogRankM(dataStartCohort, dataEndCohort, adjYELatRiskDataM, adjYELdeathsDataM, 1950)
adjYELlr1F = LogRankF(dataStartCohort, dataEndCohort, adjYELatRiskDataF, adjYELdeathsDataF, 1950)

adjTyELlr1M = LogRankM(dataStartCohort, dataEndCohort, adjTyELatRiskDataM, adjTyELdeathsDataM, 1950)
adjTyELlr1F = LogRankF(dataStartCohort, dataEndCohort, adjTyELatRiskDataF, adjTyELdeathsDataF, 1950)

paYELm = 2 * (1 - pnorm(abs(adjYELlr1M)))
paYELf = 2 * (1 - pnorm(abs(adjYELlr1F)))
paTyELm = 2 * (1 - pnorm(abs(adjTyELlr1M)))
paTyELf = 2 * (1 - pnorm(abs(adjTyELlr1F)))

# Test statistics
adjYELlr1M
adjYELlr1F
adjTyELlr1M
adjTyELlr1F

# p-values
paYELm
paYELf
paTyELm
paTyELf


# Reset data parameters
dataStartCohort = 1898
dataEndCohort = 2001



#########################################################################
#-----------------------------------------------------------------------#
#########################################################################
### Vanhuuselakeosa-calculations

# We calculate the vanhuuselakeosa using our newly estimated
# cumulative mortalities,
# and compare those to the ones given by the theory



#########################################################################
### Functions for the theorethical vanhuuselakeosa

# Constants for the interest rate,
# boundary age of the two-part model, and pension age

delta = log(1 + 0.03)
k = 70
w = 65

# D quantity for an arbitrary force of mortality

D = function(x, Lambda) {
  return(exp(-delta*x - Lambda(x)))
}

# Analytical integral of a Gompertz FoM
# In our two-part model,
# this gives the Lambda for the part specified by i (1 or 2)

LambdaIT = function(t, i, cohort, sex) {
  b = b2(cohort)
  if (i == 1) {
    if (sex == "m") {
      a1 = a1m(0, cohort)
      a2 = a2m(0, cohort)
    } else {
      a1 = a1f(0, cohort)
      a2 = a2f(0, cohort)
    }
  } else {
    if (sex == "m") {
      a1 = a1m(100, cohort)
      a2 = a2m(100, cohort)
    } else {
      a1 = a1f(100, cohort)
      a2 = a2f(100, cohort)
    }
  }
  return((a1/a2) * exp(a2*b) * (exp(a2*t) - 1))
}

# Final non-ikasiirto D quantity in our two-part model

DxTp = function(x, sex) {
  f1 = function(t) {
    return(LambdaIT(t, 1, 1950, sex))
  }
  
  f2 = function(t) {
    return(LambdaIT(t, 2, 1950, sex))
  }
  
  if (x <= k) {
    return(D(x, f1))
  } else {
    return(D(x, f2) * (D(k, f1) / D(k, f2)))
  }
}

# N quantity for a generic Lambda
# This is estimated nurerically as detailed in
# Lahti & Toro (2018), using the Simpson rule
# Basically the integral is estimated in length 1 chunks
# from x to 129 using the rule, and assumed to be 0 past 129

Nx = function(x, Lambda) {
  s = 0
  if (x %% 2 == 0) {
    for (j in seq(x, 126, by = 2)) {
      s = s + ((1/3) * (D(j, Lambda) + (4 * D(j + 1, Lambda)) + D(j + 2, Lambda)))
    }
    s = s + ((D(128, Lambda) + D(129, Lambda)) / 2)
  } else {
    for (j in seq(x, 127, by = 2)) {
      s = s + ((1/3) * (D(j, Lambda) + (4 * D(j + 1, Lambda)) + D(j + 2, Lambda)))
    }
  }
  return(s)
}

# N quantity in our specific two-part model, no ikasiirto

NxTp = function(x, sex) {
  
  f1 = function(t) {
    return(LambdaIT(t, 1, 1950, sex))
  }
  
  f2 = function(t) {
    return(LambdaIT(t, 2, 1950, sex))
  }
  
  if (x <= k) {
    return(Nx(x, f1) - Nx(k, f1) + (Nx(k, f2) * (D(k, f1) / D(k, f2))))
  } else {
    return(Nx(x, f2) * (D(k, f1) / D(k, f2)))
  }
}

# Final theorethical vanhuuselakeosa, ikasiirto included

PvT = function(x, cohort, sex) {
  b = b2(cohort)
  return(NxTp(w + b, sex) / DxTp(x + b, sex))
}



#########################################################################
### Functions for the empirical vanhuuselakeosa

# Because our data starts at 18-year olds, we need to
# account for the "lost" mortality of those first 18 years
# Also at ages above 106,
# we stop having data for some of the populations
# We use our theorethical model to correct for this,
# so in practice the first 17 Ds use theory only,
# for 18-106 we add the estimate and beyond that we use theory again

dataStartAge = 18
dataEndAge = 106

# Build the integrated FoM as described above
# We assume the NA-estimates are given in
# a vector corresponding to the ages in ageRange

getELambda = function(NAestimate, cohort, sex, ageRange) {
  s = rep(0, 129)
  for (age in 1:129) {
    currentAgeIndex = age - ageRange[1] + 1
    lastAgeIndex = dataEndAge - ageRange[1] + 1
    
    if (age < dataStartAge) {
      s[age] = LambdaIT(age, 1, cohort, sex)
    } else if (age > dataEndAge) {
      s[age] = NAestimate[lastAgeIndex] + LambdaIT(age, 2, cohort, sex) - LambdaIT(dataEndAge, 2, cohort, sex)
    } else {
      s[age] = NAestimate[currentAgeIndex] + LambdaIT(17, 1, cohort, sex)
    }
  }
  return(s)
}

# Empirical D quantity
# Only works if x is an integer

DxE = function(x, NAestimate, cohort, sex, ageRange) {
  Lambda = function(t) {
    return(getELambda(NAestimate, cohort, sex, ageRange)[t])
  }
  return(D(x, Lambda))
}

# Empirical N quantity
# Same as above, x is assumed to be an integer

NxE = function(x, NAestimate, cohort, sex, ageRange) {
  Lambda = function(t) {
    return(getELambda(NAestimate, cohort, sex, ageRange)[t])
  }
  return(Nx(x, Lambda))
}

# Empirical vanhuuselakeosa

PvE = function(x, NAestimate, cohort, sex, ageRange) {
  return(NxE(w, NAestimate, cohort, sex, ageRange) / DxE(x, NAestimate, cohort, sex, ageRange))
}



#########################################################################
### Results

ageRange = 18:w
veoYELvsTyELm = rep(0, length(ageRange))
veoYELvsTyELf = rep(0, length(ageRange))

# This one takes about 15 seconds to run
for (age in ageRange) {
  ageIndex = age - ageRange[1] + 1
  
  YELm = PvE(age, naYELm, 1950, "m", fullAgeRange)
  YELf = PvE(age, naYELf, 1950, "f", fullAgeRange)
  
  TyELm = PvE(age, naTyELm, 1950, "m", fullAgeRange)
  TyELf = PvE(age, naTyELf, 1950, "f", fullAgeRange)
  
  veoYELvsTyELm[ageIndex] = YELm / TyELm
  veoYELvsTyELf[ageIndex] = YELf / TyELf
}

# Ratios of YEL vs TyEL veo:s
plot(ageRange, veoYELvsTyELm, type = "l", xlab = "Age", ylab = "YEL/TyEL", main = "", col = "blue", log = "", ylim = c(1, 1.04))
lines(ageRange, veoYELvsTyELf, type = "l", col = "orange")
legend(x = "topright", legend=c("Miehet", "Naiset"), fill = c("blue", "orange"))



#########################################################################
# Total per capita veo:s for the two populations using
# the two different empirical mortalities

# The below takes a whopping 6 minutes for the full year range
yearRange = dataStartYear:dataEndYear
totalVeoTyEL = rep(0, length(yearRange))
totalVeoYEL1 = rep(0, length(yearRange))
totalVeoYEL2 = rep(0, length(yearRange))
sTyEL = rep(0, length(yearRange))
sYEL = rep(0, length(yearRange))

for (year in yearRange) {
  yearIndex = year - yearRange[1] + 1
  
  for (cohort in dataStartCohort:dataEndCohort) {
    age = year - cohort
    if (age >= 18 && age <= 65) {
      totalVeoTyEL[yearIndex] = totalVeoTyEL[yearIndex] + (atRisk(age, cohort, TyELatRiskDataM) * PvE(age, naTyELm, cohort, "m", fullAgeRange)) + (atRisk(age, cohort, TyELatRiskDataF) * PvE(age, naTyELf, cohort, "f", fullAgeRange))
      totalVeoYEL1[yearIndex] = totalVeoYEL1[yearIndex] + (atRisk(age, cohort, YELatRiskDataM) * PvE(age, naTyELm, cohort, "m", fullAgeRange)) + (atRisk(age, cohort, YELatRiskDataF) * PvE(age, naTyELf, cohort, "f", fullAgeRange))
      totalVeoYEL2[yearIndex] = totalVeoYEL2[yearIndex] + (atRisk(age, cohort, YELatRiskDataM) * PvE(age, naYELm, cohort, "m", fullAgeRange)) + (atRisk(age, cohort, YELatRiskDataF) * PvE(age, naYELf, cohort, "f", fullAgeRange))
      sTyEL[yearIndex] = sTyEL[yearIndex] + atRisk(age, cohort, TyELatRiskDataM) + atRisk(age, cohort, TyELatRiskDataF)
      sYEL[yearIndex] = sYEL[yearIndex] + atRisk(age, cohort, YELatRiskDataM) + atRisk(age, cohort, YELatRiskDataF)
    }
  }
  totalVeoTyEL[yearIndex] = totalVeoTyEL[yearIndex] / sTyEL[yearIndex]
  totalVeoYEL1[yearIndex] = totalVeoYEL1[yearIndex] / sYEL[yearIndex]
  totalVeoYEL2[yearIndex] = totalVeoYEL2[yearIndex] / sYEL[yearIndex]
}

dmgDiff = totalVeoYEL1/totalVeoTyEL
mrtDiff = totalVeoYEL2/totalVeoTyEL
mrtPrc = (mrtDiff - dmgDiff) / (mrtDiff - 1)

dmgDiff
mrtDiff
mrtPrc

# Plot the ratios
plot(yearRange, dmgDiff, type = "l", xlab = "Vuosi", ylab = "YEL/TyEL", main = "", col = "blue", ylim = c(1, 1.3))
lines(yearRange, mrtDiff, type = "l", col = "orange")
legend(x = "topleft", legend=c("Kuolevuus + ikarakenne", "Vain ikarakenne"), fill = c("orange", "blue"))




