# Picarro data handling attempt.  Stealing from https://github.com/bpbond/R-data-picarro/blob/master/scripts/picarro.R
library(plyr)
library(dplyr)
library(readr)     
library(ggplot2)
library(cowplot)
library(reshape2)


## Prior to running the rest of the code
## 1) Set the wd to the raw Picarro dat files

#setwd("/Users/bobhall/Dropbox/Picarro/Picarro raw data files/2018_07_13_Roys/")
###########
## IMPORT
###########


## Import meta data for each day:  In this case metadata consists of a sample identifier, and time stamp, where time stamp is when the sample syring hits about 15mL and it 
## looks like the trace is leveling off for that sample.  This metadata can and should have other stuff in it, such as sample location, collection time,
## pH, volume of sample, volume of equlibration gas (typically 70mL), temp of equilibration
meta <- read.table("2018_08_23_BlaineDiel.txt", sep=",",header=T)
meta$DATETIME <- as.POSIXct(as.character(meta$ResetTime), format="%Y-%m-%d %H:%M:%S")
meta$samplenum<-seq(1:length(meta$ResetTime))


## Import the raw data
filelist <-list.files( pattern="dat", full.names = TRUE)  ##makes a vector of that day's files
rawdata <- list()  #empty list
#below loops over all the files making a list containg a tibble for each file and stuffing it into a list
for(f in filelist) {
  cat("Reading", f, "\n") #no idea why we need this
  read_table(f) %>%
    select(DATE, TIME, ALARM_STATUS, MPVPosition, `12CO2_dry`,`13CO2`, Delta_Raw_iCO2, HP_12CH4, `HP_Delta_iCH4_Raw`) -> rawdata[[f]] #add whichever things we want here
}
rawdata <- bind_rows(rawdata)

## Fix datetime
rawdata$DATETIME <- lubridate::ymd_hms(paste(rawdata$DATE, rawdata$TIME))
rawdata$DATETIME <- as.POSIXct(as.character(rawdata$DATETIME), format="%Y-%m-%d %H:%M:%S")

# Get rid of unneeded fields
rawdata$DATE <- rawdata$TIME <- NULL

##average data for each second where we have a measure.  Not optional
rawdata<- rawdata %>% group_by(as.numeric(rawdata$DATETIME)) %>% summarise_all(funs(mean))

#Plot it. Look ok?
plot(rawdata$DATETIME,rawdata$`12CO2_dry`)
plot(rawdata$DATETIME,rawdata$Delta_Raw_iCO2)
plot(rawdata$DATETIME,rawdata$HP_12CH4)
plot(rawdata$DATETIME,rawdata$`HP_Delta_iCH4_Raw`)

##get rid of all the early in the day data to make a smaller file.  This step is optional
rawdata <- subset(rawdata, rawdata$DATETIME > "2018-09-28 20:00:00")

##if ok, then name the rawdat file for that day here
dat<-rawdata
colnames(dat)[4:7] <- c("CO2","Delta_iCO2","CH4","Delta_iCH4")

#############
## SUBSET
#############

## Subset the main datafile based on sample time + some amount of time, 
## 40 seconds appears correct, but always check this time.  
## To short and we lost data.  Too long and we average off the plateau
meta$ENDTIME<-meta$DATETIME+30

subdat<-list()
for(i in 1:length(meta$samplenum)){
  
subdat[[i]]<-subset(dat, DATETIME>meta$DATETIME[i]  & DATETIME< meta$ENDTIME[i])
subdat[[i]]$samplenum <-i
}
subdat <- bind_rows(subdat)

## Plot
subdat_melt <- melt(subdat[,4:9], id.vars=c("DATETIME","samplenum"))
ggplot(subdat_melt, aes(DATETIME, value, color=samplenum))+
  geom_point()+
  facet_wrap(~variable, ncol=2, scales = "free")


## Summarize per sample
subdatsum<- subdat %>% group_by(samplenum) %>% summarise(CO2=mean(CO2), delCO2=mean(Delta_iCO2), 
                                                         CH4=mean(CH4), delCH4=mean(Delta_iCH4) )


##merge with meta, and give it a new name

#Thinking about the metafile, here is what should also go in
#water_temp, equil_temp, pH, ANC, air_vol (it won't always be 70), baro_press (absolute, mm Hg)
#not the  kH_CO2 because that is a derived value that we will calc each time

sum_file<- full_join(meta,subdatsum)  #there you go, joining by samplenum


## Remove standards
sum_file <- na.omit(sum_file)

## Remove any problematic syringe lines based on notes
#sum_file <- sum_file[-which(sum_file$SampleID == "Blaine_NoClearSignal"),]
sum_file <- sum_file[-which(sum_file$Syringe %in% c(28,95)),]


## Plot
sum_file_melt <- melt(sum_file[,c("CO2","delCO2","CH4","delCH4","DATETIME","samplenum")], id.vars=c("DATETIME","samplenum"))
ggplot(sum_file_melt, aes(DATETIME, value, color=samplenum))+
  geom_point()+
  facet_wrap(~variable, ncol=2, scales = "free")



## merge again with other data if available
notes <- read.table("2018_09_28_OakCreek_notes.txt", sep=",",header=T)
notes <- notes %>% 
  separate(Syringes, c("SyringeA","SyringeB","SyringeC"))


sumA <- merge(sum_file, notes, by.x="Syringe", by.y = "SyringeA")
sumB <- merge(sum_file, notes, by.x="Syringe", by.y = "SyringeB")
sumC <- merge(sum_file, notes, by.x="Syringe", by.y = "SyringeC")
sum_file <- rbind(sumA[,1:17], sumB[,1:17], sumC[,1:17])
## Fix names
sum_file <- sum_file[,c(1:4,9:17)]
colnames(sum_file)[9] <- "SampleID"

## OR
#sum_file <- merge(sum_file, notes, by="SampleID")


############
## Export
###########
write.csv(sum_file, "2018_09_28_OakDielprocessed.csv")



#########################
## CALCULATE CO2 IN H20 (in the post-processing files now)
#########################
#I (Bob) think this is correct

CO2calc2<-function(x) {
  x$kH_CO2<-29* exp( 2400 * ((1/(x$equiltemp+273.15)) - 1/(298.15)) ) #temp correction for Henry's law
  x$CO2_mol_water<- x$CO2_vol_air  / x$kH_CO2
  x$CO2_mol_air<- x$CO2_vol_air/ ((x$equiltemp+273.15)*0.08026/x$pressure)
  
  x$CO2_conc<-(x$CO2_mol_water*x$volwater + x$CO2_mol_air*x$volair)/100 # µmol/L or mmol/m3
  return(x)
  
}


## Need to set certain columns to desired values if all the same
## If each sample is different, simple set the column to a column specified in the metadata file
names(sum_file)

sum_file$CO2_vol_air <- sum_file$CO2
sum_file$equiltemp <- sum_file$WaterTemp_C
sum_file$volair <- sum_file$air_mL
sum_file$volwater <- sum_file$H2O_mL
sum_file$pressure <- 0.9

## Run function on datafile
sum_file <- CO2calc2(sum_file)
names(sum_file)

ggplot(sum_file, aes(SampleID, CO2_conc))+
  geom_point(size=3)+
  theme(axis.text.x=element_text(angle=45, hjust = 1))





############
## OLD
#############
###Warning, the below is incorrect

# # Calculate CO2 -- see details below
# CO2calc<-function(x) {
#   x$kH_CO2<-0.034* exp( 1300 * ((1/(x$equiltemp+273.15)) - 1/(298.15)) ) #temp correction for Henry's law
#   x$Hcc_CO2<- x$kH_CO2 *0.08205*(x$equiltemp+273.15)
#   x$CO2_vol_water <-  x$Hcc_CO2 * x$CO2_vol_air
#   x$CO2_mol_water<- x$CO2_vol_water / ((x$equiltemp+273.15)*0.08026/x$pressure)
#   x$CO2_mol_air<- x$CO2_vol_air/ ((x$equiltemp+273.15)*0.08026/x$pressure)
#   
#   x$CO2_conc<-x$CO2_mol_water + (x$CO2_mol_air*x$volair/x$volwater) # µmol/L or mmol/m3
#   return(x)
#   
# }





##Then use Henry's law to calculate the concentration of CO2 in the water based on the conc in the headspace. Hcc_CO2 is the ratio of conc in the water to conc in the air and 
#is therefore unitless.  tempK is temp in K (273.15+degC)

## If not already in meta file and approximately the same, specify the temperature of all samples
subdat$tempK <- 273.15+21

subdat$kH_CO2<-0.034* exp( 1300 * ((1/subdat$tempK) - 1/(298.15)) ) #temp correction for Henry's law
subdat$Hcc_CO2<- subdat$kH_CO2 *0.08205*subdat$tempK #convert pressure based henry's law to conc based.


## For water, 
subdat$CO2_vol_air <- subdat$CO2

subdat$CO2_vol_water <-  subdat$Hcc_CO2 * subdat$CO2_vol_air #(the latter term comes from the machine)
#Need to convert concnetration that the machine measures (ppmv = µL/L) to µmol/L.  Rememeber that PV=nRT and R= 0.08206 L*atm/(K*mol). 
#P at station = 0.9 atm
subdat$CO2_mol_water<- (subdat$tempK*0.08026/0.9)*subdat$CO2_vol_water 
subdat$CO2_mol_air<- (subdat$tempK*0.08026/0.9)*subdat$CO2_vol_air

##now calculate the total CO2 by weighting for water and air volumes and summing
## If not already specified in meta file, specify volair and volwater
subdat$volair <- 0.07
subdat$volwater <- subdat$H2O_mL

subdat$CO2totwat<-subdat$CO2_mol_water + (subdat$CO2_mol_air*subdat$volair/subdat$volwater)






