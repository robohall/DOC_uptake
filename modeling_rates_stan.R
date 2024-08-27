##
#modeling rates, but with Stan and not one by one

library (ggplot2)
library(rstan)
library(tidybayes)
#get processed data

data_doc<- read.csv("./data_doc.csv")
site_info<-read.csv("./site_data.csv")
source("./doc_functions.R")

data_doc$trial<- as.integer(factor(data_doc$code, levels = c("b_g1_up","b_g1_down","b_g2_up","b_g2_down","b_l_up","b_l_down")))


KQ_out<- function(data, site_data){
  Q<-site_data$salt_mass*2100/nutint(data$MinFrom0, data$SpC_corr, bkg=0) 
  int13c<- nutint(data$MinFrom0, data$doc_13C, bkg=0 ) 
  int13c_flux<-int13c*Q
  K<-  -(log(int13c_flux) -log(site_data$doc_mass_13C_mg))/site_data$length
  velocity<- site_data$length/site_data$travel_time
  z<-0.001*Q/(site_data$width*velocity)
  

  
  c(K, Q, velocity, z)
  
}



KQ<-data.frame(K=NA, q=NA, velocity=NA, z=NA)

for ( i in 1:6){
  
  dat<- data_doc[data_doc$trial==i, ]
  KQ[i,]<- KQ_out(data=dat, site_data=site_info[i,] )
  
}

code<-c("b_g1_up","b_g1_down","b_g2_up","b_g2_down","b_l_up","b_l_down")
KQ$code <-code

KQ$vf<-( KQ$q*0.001*1440/site_info$width[1:6])* KQ$K

vf_mmmin <- KQ$vf*1000/1440 ##mm/min
mean(vf_mmmin[1:4])

##save lots of stan coding by --calculating--R_g outside of it
R_g_out<- function(data,Q){
  removed_c_flux = data$removed_c*Q
  
  r_g<-0
  for (i in 2: length (data$MinFrom0) ){
    int_u<- (data$MinFrom0[i]-data$MinFrom0[i-1])*(removed_c_flux[i])
    r_g[i]<- r_g[i-1]+int_u
    
  }
  r_g
}

R_g<-NA
for ( i in 1:6){
  
  dat<- data_doc[data_doc$trial==i, ]
 R_g_one<- R_g_out(data=dat, Q=KQ$q[i])
 R_g<-c(R_g,R_g_one)
  
}
R_g<-R_g[-1]

data_doc$R_g<-R_g

data_doc <- left_join(data_doc,KQ, by="code")


##data for stan model

dic_dat<- list(N=length(data_doc$trial), N_groups=6, DIC_flux = data_doc$DIC13 * data_doc$q, R_g=data_doc$R_g,
               MinFrom0 = data_doc$MinFrom0, group_id = data_doc$trial, removed_c = data_doc$removed_c, q = data_doc$q)


dic_fit<-stan(file= "dic_prod.stan", data=dic_dat,  iter = 2000, chains = 4)


##extract params and samples to make a figure

doc_steps<-  dic_fit %>% 
  spread_draws(kc[i], fi[i])

med_doc_steps<- median_qi(doc_steps)

head(doc_steps)
  
pred_steps<- dic_fit %>% 
  spread_draws(pred_imm[i], pred_delayed[i], ndraws=200)
pred_steps$group<- rep(data_doc$trial, each=200)
pred_steps$MinFrom0<- rep(data_doc$MinFrom0, each=200)
pred_steps$MinFrom0<- rep(data_doc$MinFrom0, each=200)
pred_steps$step<-rep(1:217, 200)
pred_steps$pred_tot <- pred_steps$pred_imm + pred_steps$pred_delayed

pred_steps <-pred_steps %>% arrange(step, i)

head(pred_steps)


###extract

####
#####
#Some figures, final drafts are in plots file, delete these at some point

library(wesanderson)

ggplot( doc_steps, aes(x=kc, y=fi, color = as.factor(i))) +
  geom_point(alpha=0.5, size =0.1) + ylim(0,0.8) +xlim(0,8)+
  theme_classic()+
  scale_colour_manual(values=c("#EBCC2A","#F2AD00","#D69C4E", "#F98400", "#ABDDDE","#5BBCD6"),name="", breaks=c(1,2,3,4,5,6),
                      labels=c("g_1_up", "g_1_down", "g_2_up", "g_2_down", "l_up", "l_down"))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  labs(x = expression(paste(italic(k[c])," (1/d)")), y = expression(italic(f[i])))

##1/kc

ggplot( doc_steps) +
  geom_point( aes(x=1/kc, y=fi, color = as.factor(i)),alpha=0.1, size =0.1) + ylim(0,0.8) +xlim(0,1)+
  theme_classic()+
  scale_colour_manual(values=c("#EBCC2A","#F2AD00","#D69C4E", "#F98400", "#ABDDDE","#5BBCD6"),name="", breaks=c(1,2,3,4,5,6),
                      labels=c("g_1_up", "g_1_down", "g_2_up", "g_2_down", "l_up", "l_down"))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  labs(x = expression(paste(italic(1/k[c])," (d)")), y = expression(italic(f[i])))+
  geom_point(data=med_doc_steps, aes(x=1/kc, y=fi))



##Uptake length and Vf summaries fo paper
KQ$sw<- 1/ KQ$K
KQ$vf<- KQ$K * KQ$z *KQ$velocity *1440  # units m/d

mean(KQ$sw[1:4])
mean(KQ$sw[5:6])

mean(KQ$vf[1:4])
mean(KQ$vf[5:6])
