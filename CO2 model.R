## DIC forward model

# $[CO_{2}]$ = C
# $[HCO_{3}^-]$ = B
# $[CO_{3}^{2-}]$ = Ca
# Alkalinity = A
# DIC = D

##Date, trimmed

pco2m<-pco2[109:972,]
dat<-pco2m

alk<-5.3

################################
## Functions
###############################

## Equilibrium constants
K1calc<- function(temp) { 10^( (-3404.71/(273.15+temp)) + 14.844 -0.033*(temp+273.15) )}
K2calc<- function(temp) { 10^( (-2902.39/(273.15+temp)) + 6.498 -0.0238*(temp+273.15) )}

##need CO2 saturation equation
Khcalc<- function(temp){
  A1<- -58.0931
  A2<- 90.5069
  A3<- 22.2940
  exp( A1 + A2*100/(273.15+temp) + A3*log( (273.15+temp)/100) )
  
}
Khcalc(25)
1/Khcalc(25)


#Calculate CO2 sat in mmol/l
csat<- function(temp,bp) {
  c<- 400e-3*(bp/760)*Khcalc(temp)
  c
}
csat(25,760)

## Carbonate Chemistry from C & A
Carbfrom_C_A <- function(K1, K2, C, A){
  H_minus <- (((-K1*C))-sqrt(((K1*C)^2)-(4*-1*A*2*K1*K2*C)))/(2*-1*A)
  pH <- -1*log10(H_minus)
  B <- (K1*C)/H_minus
  Ca <- (K2*B)/H_minus
  D <- C + B + Ca
  A_check <- B + 2*Ca
  
  l <- list(H_minus, pH, C, B, Ca, D, A, A_check)
  names(l) <- c("H", "pH", "C", "B", "Ca", "D", "A","A check")
  return(l)
}


## Carbonate Chemistry from D & A


Carbfrom_D_A <- function(K1, K2, D1, A){
  a <- A
  b <- K1*(A-D1)
  c <- (A-(2*D1))*K1*K2
  H_t <- ((-1*b)+sqrt((b^2)-(4*a*c)))/(2*a)
  
  pH_t <- -1*log10(H_t)
  
  B_t <- (D1*K1*H_t)/((H_t^2)+(K1*H_t)+(K1*K2))
  Ca_t <- (D1*K1*K2)/((H_t^2)+(K1*H_t)+(K1*K2))
  
  C_t <- (H_t*B_t)/K1
  D2 <- C_t + B_t + Ca_t
  A2_check <- B_t + 2*Ca_t
  
  l <- list(H_t, pH_t, C_t, B_t, Ca_t, D2, A, A2_check)
  names(l) <- c("H", "pH", "C", "B", "Ca", "D", "A", "A check")
  return(l)
}


Carbfrom_D_A(K1, K2, 5.473, A1_molL)


## This function calculates KO2 at any given tempertaure from K600. via schmidt number scaling.  The scaling equation if From Jähne et al. (1987), Schmidt number conversions from Wanninkhof et al. 1992.
Kcor<-function (temp,K600) {
  K600/(600/(1800.6-(temp*120.1)+(3.7818*temp^2)-(0.047608*temp^3)))^-0.5
}

###########################################################################
## Start by establishing the first DIC value using CO2 and alkalinity
###########################################################################
# Specify SampleID in pH_Alk_data


# Need K1, K2, CO2 mol/L, and Alkalinity mol/L for Carbfrom_C_A function

# Calc starting K1 and K2
K1 <- K1calc(dat$white_temp[1])
K2 <- K2calc(dat$white_temp[1])

# Calc starting A in mol/L
#A1 <- pH_Alk_data[which(pH_Alk_data$SampleID == start_ID),]$Alk_mgLCaCO3
#A1_molL <- A1*(1/100.0869)*(1/1000) ## actually convert to mol/L
A1_molL<-alk
D1_CarbEq <- Carbfrom_C_A(K1, K2, dat$white_co2_conc[1]/12, A1_molL)
D1_CarbEq

##########################################################################
## Next calculate the next DIC time step
##########################################################################
## Step 2: Create a forward moving model of DIC

#D2 = D1 + ((GPP/Z)\times(PAR/sum(PAR)) + (ER/z)\times dt + K(Csat - C1) + GW\times dt

## Where:
# GPP, ER, and K are from the oxygen model (mol/L)
# Temperature and light are directly from data
# Csat is from Henry's law
# C1 is from the previous timestep
# GW will be calculated as an anomaly

GPP <- -10/32 ##  mol m2 d
ER <- 11/32
K <- 13.5
z <- 0.2
ts <- 10/(60*24)
bp <- 0.9*760
light <- dat$neplight
temp <- dat$white_temp


D.mod <- numeric(length(dat$white_co2))
D.mod[1]<-D1_CarbEq$D
C.mod<-dat$white_co2_conc[1]/12
# forward model

##this is wrong!


for (i in 2:length(dat$white_temp)) {
  
  D.mod[i]<- D.mod[i-1] + 
    ((GPP/z)*(light[i]/(sum(light)/6))) + 
    ER*ts/z +
    (Kcor(temp[i],K))*ts*(csat(temp[i],bp) - C.mod[i-1]) 
Ceq<-Carbfrom_D_A(D1=D.mod[i-1],K1=K1calc(temp[i-1]), K2=K2calc(temp[i-1]), A=alk)
C.mod[i]<-Ceq$C
}

plot(C.mod*12)

pco2m$C.mod<- C.mod*12

pco2m$C.mod.flux<-KCO2fromK600(pco2m$white_temp,k_blaine) *(pco2m$C.mod-pco2m$satconc)
plot(pco2m$C.mod.flux)

# visualize
plot(seq(1:length(dat$CO2_mmolL)),
     D.mod, type="l",xlab="Time", ylab="CO2 (mg/L)",
     cex.lab=1.5, cex.axis=1.5, lwd=2 )
points(seq(1:length(dat$CO2_mmolL)), dat$CO2_mmolL)

