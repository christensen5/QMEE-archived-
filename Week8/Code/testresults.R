setwd("/home/alexander/Documents/QMEE/Week8/Results/hpcrun")
filenames = list.files(getwd(), pattern = "AKC_NTSresult_iter*")
numResults = length(filenames)

for (i in seq(1,numResults)) {
  loadstr = paste0(getwd() , "/" , filenames[i])
  load(loadstr)
  print(paste0(filenames[i],", ",length(abundances)))
}
