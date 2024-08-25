library(tidyverse)
library(lubridate)

setwd("/Users/phill/Dropbox/Isotopes/sugarIsotopes/201204")

dat <- read_csv("13CAdditionData_MovieWellGlucose_20190810_analyzed.csv")
picarro <- read_csv("CFIDS2157_IsotopicData_20210119_224851.csv")
aurora <- read_csv("Hall_210119.csv")
doc <- read_csv("doc.csv")

inds <- c() #indices of samples in original data file

for(SampleID in doc$SampleID){
  inds <- c(inds, match(SampleID, dat$DOC.ID))
}
inds <- na.omit(inds)

set<-dat[inds,] #Pick out rows of the samples we ran from the original data file

#Attach DOC values from the DOC file
for (i in 1:length(set$DOC.ID)){
  set$DOC_conc[i] <- doc$DOC[match(set$DOC.ID[i], doc$SampleID)]
}


auroraSampleNumbers <- c() #Corresponding "Spl#" in the Aurora output of the samples in the current dat file

for(SampleID in set$DOC.ID){
  auroraSampleNumbers <- c(auroraSampleNumbers, aurora$`Spl#`[match(SampleID, aurora$`Sample ID`)])
}

samplePicarroIndices <- c() #Indices of the samples we ran in the Picarro output
for (n in auroraSampleNumbers) {
  samplePicarroIndices <-
    c(samplePicarroIndices, match(n, picarro$`Sample Id`))
}

#Make list of Delta values of the samples from the Picarro output
deltas <- c()
for (i in samplePicarroIndices) {
  deltas <- c(deltas, picarro$`Delta 13C`[i])
}

set$`Delta 13C` <- deltas #Attach delta values
set$Datetime <-
  mdy_hms(paste(as.character(set$Date), as.character(set$SampleTime))) #Attach datetime for plots

plot(set$Datetime, set$`Delta 13C`)
plot(set$Datetime, set$DOC_conc)
plot(1:length(set$DOC_conc), set$DOC_conc)
plot(set$DOC_conc, set$`Delta 13C`)

#Write DOC and Deltas back to original file
for(i in 1:length(set$DOC.ID)){
  dat$DOC_conc[match(set$DOC.ID[i],dat$DOC.ID)] <- set$DOC_conc[i]
  dat$`Delta 13C`[match(set$DOC.ID[i],dat$DOC.ID)] <- set$`Delta 13C`[i]
}

write.csv(dat, file="13CAdditionData_MovieWellGlucose_20190810_analyzed.csv")
