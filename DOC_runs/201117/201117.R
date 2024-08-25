library(tidyverse)
library(lubridate)

setwd("/Users/phill/Dropbox/Isotopes/sugarIsotopes/201117")

dat <- read_csv("13CAdditionData_CrestonLeachate_20190809(1).csv") #original data file from experiment
picarro <- read_csv("CFIDS2157_IsotopicData_20201117_183150.csv") #Picarro output
aurora <- read_csv("Hall_201117.csv") #Aurora output
doc <- read_csv("doc.csv") #Calculate DOC from Aurora output

inds <- c() #indices of samples in original data file
auroraSampleNumbers <- c() #Correspondings "Spl#" of the samples in the Aurora output
for(SampleID in doc$SampleID){
  inds <- c(inds, match(SampleID, dat$DOC.ID))
  auroraSampleNumbers <- c(auroraSampleNumbers, aurora$`Spl#`[match(SampleID, aurora$`Sample ID`)])
}

set<-dat[inds,] #Pick out rows of the samples we ran from the original data file
set$DOC <- doc$DOC #Attach DOC values from the DOC file

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

set$Delta <- deltas #Attach delta values
set$Datetime <-
  mdy_hms(paste(as.character(set$SampleDate), as.character(set$SampleTime))) #Attach datetime for plots

plot(set$Datetime, set$Delta)
plot(set$Datetime, set$DOC)
plot(1:length(set$DOC), set$DOC)
plot(set$DOC, set$Delta)

dat$DOC <- NA
dat$Delta <- NA
#Write DOC and Deltas back to original file
for(i in 1:length(set$DOC.ID)){
  dat$DOC[match(set$DOC.ID[i],dat$DOC.ID)] <- set$DOC[i]
  dat$Delta[match(set$DOC.ID[i],dat$DOC.ID)] <- set$Delta[i]
}

write.csv(dat, file="13CAdditionData_CrestonLeachate_20190809(1)_Analyzed.csv")
