## Load mean temperature data for Key West, Florida, and perform statistical correlation analysis.

load('../Data/KeyWestAnnualMeanTemperature.RData') #load data

datamat = as.matrix(ats) #store at matrix for efficiency

natCor = cor(datamat[1:99,2], datamat[2:100,2]) #compute correlation between pairs of successive years

numPerms = 10000
permCors = matrix(NA,numPerms,1) #pre-allocate matrix
for (i in 1:numPerms) {
  perm = sample(datamat[,2]) #permute years randomly
  permCors[i] = cor(perm[1:99], perm[2:100]) #compute correlation between pairs of successive (post-permutation) years
}

p.eff = sum(permCors > natCor)/numPerms #effective P-value is the proportion of post-permutation correlation coeffs which are greater than initial correlation coeff

print(paste("The effective P-value is", as.character(p.eff)))