####  
library(rethinking)
library(data.table)
rm(list=ls(all=TRUE))




data<-fread("./data/AllData.csv", header=T, colClasses="character", na.strings = "")
names(data)

  data$TOPIC = factor(data$TOPIC)
levels(data$TOPIC)


levels(data$TOPIC) <- c("AD", "CT", "DD", "IC", "MI", "PC", "PP", "SS", "SP")

data$TAXON = factor(data$TAXON)
levels(data$TAXON)
data$TAXON[which(data$TOPIC == "PP")] = "O"


#get the number of publications per topic, region, taxon

N.data = data[, .N, by=c("TOPIC", "REGION", "TAXON")]


# get the model matrix structure
(mod.mat <- as.data.frame(model.matrix(~ REGION*TOPIC, N.data)))


# get the bayesian model 
glimmer(N ~ REGION*TOPIC + (1|TAXON), family = poisson, N.data)

#Change the names for easier reading
names(mod.mat)[2:10] <-c("Trop", "CT", "DD", "IC", "MI", "PC", "PP", "SS", "SP")
names(mod.mat)[11:18] <- c("Trop_X_CT", "Trop_X_DD", "Trop_X_IC", "Trop_X_MI", "Trop_X_PC",
  "Trop_X_PP", "Trop_X_SS", "Trop_X_SP")


# make the dataset for the model
d <-as.data.frame(cbind(N.data[,c("N","TAXON") ], mod.mat[,2:18]))
names(d)[1] <- "Npapers"
levels(d$TAXON)<- 1: length(levels(d$TAXON))
d$TAXON <- as.numeric(d$TAXON)


# model with interactions

model.1 <- map2stan(
  
  alist(
    Npapers ~ dpois( lambda ),
    log(lambda) <- Intercept +
      b_Trop*Trop +
      b_CT*CT +
      b_DD*DD +
      b_IC*IC +
      b_MI*MI +
      b_PC*PC +
      b_PP*PP +
      b_SS*SS +
      b_SP*SP +
      b_Trop_X_CT*Trop_X_CT +
      b_Trop_X_DD*Trop_X_DD +
      b_Trop_X_IC*Trop_X_IC +
      b_Trop_X_MI*Trop_X_MI +
      b_Trop_X_PC*Trop_X_PC +
      b_Trop_X_PP*Trop_X_PP +
      b_Trop_X_SS*Trop_X_SS +
      b_Trop_X_SP*Trop_X_SP +
      v_Intercept[TAXON],
    Intercept ~ dnorm(0,10),
    b_Trop ~ dnorm(0,10),
    b_CT ~ dnorm(0,10),
    b_DD ~ dnorm(0,10),
    b_IC ~ dnorm(0,10),
    b_MI ~ dnorm(0,10),
    b_PC ~ dnorm(0,10),
    b_PP ~ dnorm(0,10),
    b_SS ~ dnorm(0,10),
    b_SP ~ dnorm(0,10),
    b_Trop_X_CT ~ dnorm(0,10),
    b_Trop_X_DD ~ dnorm(0,10),
    b_Trop_X_IC ~ dnorm(0,10),
    b_Trop_X_MI ~ dnorm(0,10),
    b_Trop_X_PC ~ dnorm(0,10),
    b_Trop_X_PP ~ dnorm(0,10),
    b_Trop_X_SS ~ dnorm(0,10),
    b_Trop_X_SP ~ dnorm(0,10),
    v_Intercept[TAXON] ~ dnorm(0,sigma_TAXON),
    sigma_TAXON ~ dcauchy(0,1)
  ), data = d, chains = 4, iter = 5000
  
  
)

#Check if the model converges
tracerplot(model.1)

# Model without interactions
model.2 <- map2stan(
  
  alist(
    Npapers ~ dpois( lambda ),
    log(lambda) <- Intercept +
      b_Trop*Trop +
      b_CT*CT +
      b_DD*DD +
      b_IC*IC +
      b_MI*MI +
      b_PC*PC +
      b_PP*PP +
      b_SS*SS +
      b_SP*SP +
      v_Intercept[TAXON],
    Intercept ~ dnorm(0,10),
    b_Trop ~ dnorm(0,10),
    b_CT ~ dnorm(0,10),
    b_DD ~ dnorm(0,10),
    b_IC ~ dnorm(0,10),
    b_MI ~ dnorm(0,10),
    b_PC ~ dnorm(0,10),
    b_PP ~ dnorm(0,10),
    b_SS ~ dnorm(0,10),
    b_SP ~ dnorm(0,10),
    v_Intercept[TAXON] ~ dnorm(0,sigma_TAXON),
    sigma_TAXON ~ dcauchy(0,1)
  ), data = d, chains = 4, iter = 5000, WAIC = T
  
  
)

#check if the model converges
tracerplot(model.2)

# compare the models, is the interaction term important
compare(model.1, model.2)



# get the parameters
(sumStat <- precis(model.1, depth = 2, prob = .95, digits = 5))
sumStat[sumStat[,3]>0,]
sumStat[sumStat[,4]<0,]

write.csv(sumStat, "Parameters_Npapers.csv")


# Get the posterior probabilities and plot the values for each treatment combination
pos <- extract.samples(model.1)

graphics.off()



 mu.link <- function(Trop, CT, DD, IC, MI, PC, PP, SS, SP, 
                     Trop_X_CT, Trop_X_DD, Trop_X_IC, Trop_X_MI, Trop_X_PC, Trop_X_PP, Trop_X_SS, Trop_X_SP
                     
                     ){

   Y.est <- with(pos, Intercept + b_Trop*Trop + b_CT*CT + b_DD*DD + b_IC*IC + 
                   b_MI*MI + b_PC*PC + b_PP*PP + b_SS*SS+ b_SP*SP + 
                   b_Trop_X_CT*Trop_X_CT + b_Trop_X_DD*Trop_X_DD + b_Trop_X_IC*Trop_X_IC + 
                   b_Trop_X_MI*Trop_X_MI + b_Trop_X_PC*Trop_X_PC + b_Trop_X_PP*Trop_X_PP + 
                   b_Trop_X_SS*Trop_X_SS+ b_Trop_X_SP*Trop_X_SP 
                   )
   return(exp(Y.est))

 }

 

 op<-par(mfrow=c(3,3),mar=c(3.5,1.5,1.5,2))
 
 
# Adaptation x region  
 mu.Trop_X_AD <- mu.link(Trop = 1, CT = 0,DD = 0,IC = 0,MI = 0,PC = 0,PP = 0,SS = 0,SP = 0, 
                         Trop_X_CT = 0,Trop_X_DD = 0,Trop_X_IC = 0,Trop_X_MI = 0,Trop_X_PC = 0,
                         Trop_X_PP = 0, Trop_X_SS = 0, Trop_X_SP = 0
                         )

 mu.Tem_X_AD <- mu.link(Trop = 0, CT = 0,DD = 0,IC = 0,MI = 0,PC = 0,PP = 0,SS = 0,SP = 0, 
                         Trop_X_CT = 0,Trop_X_DD = 0,Trop_X_IC = 0,Trop_X_MI = 0,Trop_X_PC = 0,
                         Trop_X_PP = 0, Trop_X_SS = 0, Trop_X_SP = 0
 )

 dens( mu.Trop_X_AD, show.HPDI = .95, col="red",  xlim= c(0,30 ), xlab = "N papers", main = "Adaptation")
 dens( mu.Tem_X_AD, show.HPDI = .95, add = T )


# Climate Tolerance x region  

mod.mat[which(mod.mat$Trop == 1 & mod.mat$CT ==1),]

mu.Trop_X_CT <- mu.link(Trop = 1,        CT = 1,       DD = 0,              IC = 0,       MI = 0,      PC = 0,        PP = 0,        SS = 0,       SP = 0, 
                                  Trop_X_CT = 1,Trop_X_DD = 0,Trop_X_IC = 0,Trop_X_MI = 0,Trop_X_PC = 0,Trop_X_PP = 0, Trop_X_SS = 0, Trop_X_SP = 0
)

mod.mat[which(mod.mat$Trop == 0 & mod.mat$CT ==1),]
mu.Tem_X_CT <- mu.link(Trop = 0,        CT = 1,       DD = 0,              IC = 0,       MI = 0,      PC = 0,        PP = 0,        SS = 0,       SP = 0, 
                                Trop_X_CT = 0,Trop_X_DD = 0,Trop_X_IC = 0,Trop_X_MI = 0,Trop_X_PC = 0,Trop_X_PP = 0, Trop_X_SS = 0, Trop_X_SP = 0
)

dens( mu.Trop_X_CT, show.HPDI = .95, col="red", xlab = "N papers", main = "Climate tolerance", xlim=c(0,20))
dens( mu.Tem_X_CT, show.HPDI = .95, add = T) 


# density x region  

mod.mat[which(mod.mat$Trop == 1 & mod.mat$DD ==1),]

mu.Trop_X_DD <- mu.link(Trop = 1,        CT = 0,       DD = 1,              IC = 0,       MI = 0,      PC = 0,        PP = 0,        SS = 0,       SP = 0, 
                        Trop_X_CT = 1,Trop_X_DD = 1,Trop_X_IC = 0,Trop_X_MI = 0,Trop_X_PC = 0,Trop_X_PP = 0, Trop_X_SS = 0, Trop_X_SP = 0
)

mod.mat[which(mod.mat$Trop == 0 & mod.mat$DD ==1),]
mu.Tem_X_DD <- mu.link(Trop = 0,        CT = 0,       DD = 1,              IC = 0,       MI = 0,      PC = 0,        PP = 0,        SS = 0,       SP = 0, 
                       Trop_X_CT = 0,Trop_X_DD = 0,Trop_X_IC = 0,Trop_X_MI = 0,Trop_X_PC = 0,Trop_X_PP = 0, Trop_X_SS = 0, Trop_X_SP = 0
)

dens( mu.Tem_X_DD, show.HPDI = .95, xlab = "N papers", main = "Density dependence", xlim=c(0,20)) 
dens( mu.Trop_X_DD, show.HPDI = .95, col="red", add = T)



# Interspecific comp x region  

mod.mat[which(mod.mat$Trop == 1 & mod.mat$IC ==1),]

mu.Trop_X_IC <- mu.link(Trop = 1,        CT = 0,       DD = 0,              IC = 1,       MI = 0,      PC = 0,        PP = 0,        SS = 0,       SP = 0, 
                        Trop_X_CT = 1,Trop_X_DD = 1,Trop_X_IC = 1, Trop_X_MI = 0,Trop_X_PC = 0,Trop_X_PP = 0, Trop_X_SS = 0, Trop_X_SP = 0
)

mod.mat[which(mod.mat$Trop == 0 & mod.mat$IC ==1),]
mu.Tem_X_IC <- mu.link(Trop = 0,        CT = 0,       DD = 0,              IC = 1,       MI = 0,      PC = 0,        PP = 0,        SS = 0,       SP = 0, 
                       Trop_X_CT = 0,Trop_X_DD = 0,Trop_X_IC = 0,Trop_X_MI = 0,Trop_X_PC = 0,Trop_X_PP = 0, Trop_X_SS = 0, Trop_X_SP = 0
)

dens( mu.Tem_X_IC, show.HPDI = .95, xlab = "N papers", main = "Interspecific competition", xlim=c(0,50)) 
dens( mu.Trop_X_IC, show.HPDI = .95, col="red",add = T)



# Mimicry x region  


mu.Trop_X_MI <- mu.link(Trop = 1,        CT = 0,       DD = 0,              IC = 0,       MI = 1,      PC = 0,        PP = 0,        SS = 0,       SP = 0, 
                        Trop_X_CT = 0,Trop_X_DD = 0,Trop_X_IC = 0, Trop_X_MI = 1,Trop_X_PC = 0,Trop_X_PP = 0, Trop_X_SS = 0, Trop_X_SP = 0
)


mu.Tem_X_MI <- mu.link(Trop = 0,        CT = 0,       DD = 0,              IC = 0,       MI = 1,      PC = 0,        PP = 0,        SS = 0,       SP = 0, 
                       Trop_X_CT = 0,Trop_X_DD = 0,Trop_X_IC = 0,Trop_X_MI = 0,Trop_X_PC = 0,Trop_X_PP = 0, Trop_X_SS = 0, Trop_X_SP = 0
)

dens( mu.Tem_X_MI, show.HPDI = .95, xlab = "N papers", main = "Mimicry", xlim=c(0,20)) 
dens( mu.Trop_X_MI, show.HPDI = .95, col="red",add = T)



### PARENTAL CARE
mu.Trop_X_PC <- mu.link(Trop = 1,        CT = 0,       DD = 0,              IC = 0,       MI = 0,      PC = 1,        PP = 0,        SS = 0,       SP = 0, 
                        Trop_X_CT = 0,Trop_X_DD = 0,Trop_X_IC = 0, Trop_X_MI = 0,Trop_X_PC = 1,Trop_X_PP = 0, Trop_X_SS = 0, Trop_X_SP = 0
)


mu.Tem_X_PC <- mu.link(Trop = 0,        CT = 0,       DD = 0,              IC = 0,       MI = 0,      PC = 1,        PP = 0,        SS = 0,       SP = 0, 
                       Trop_X_CT = 0,Trop_X_DD = 0,Trop_X_IC = 0,Trop_X_MI = 0,Trop_X_PC = 0,Trop_X_PP = 0, Trop_X_SS = 0, Trop_X_SP = 0
)
dens( mu.Trop_X_PC, show.HPDI = .95, col="red", xlab = "N papers", main = "Parental care", xlim=c(1,30)) 
dens( mu.Tem_X_PC, show.HPDI = .95,add = T)


### Predator prey 
mu.Trop_X_PP <- mu.link(Trop = 1,        CT = 0,       DD = 0,              IC = 0,       MI = 0,      PC = 0,        PP = 1,        SS = 0,       SP = 0, 
                        Trop_X_CT = 0,Trop_X_DD = 0,Trop_X_IC = 0, Trop_X_MI = 0,Trop_X_PC = 0,Trop_X_PP = 1, Trop_X_SS = 0, Trop_X_SP = 0
)


mu.Tem_X_PP <- mu.link(Trop = 0,        CT = 0,       DD = 0,              IC = 0,       MI = 0,      PC = 1,        PP = 0,        SS = 0,       SP = 0, 
                       Trop_X_CT = 0,Trop_X_DD = 0,Trop_X_IC = 0,Trop_X_MI = 0,Trop_X_PC = 0,Trop_X_PP = 0, Trop_X_SS = 0, Trop_X_SP = 0
)
dens( mu.Tem_X_PP, show.HPDI = .95, xlab = "N papers", main = "Predator-prey interactions", xlim =c(1,80)) 
dens( mu.Trop_X_PP, show.HPDI = .95, col="red",add = T)


### Sexual selection 
mu.Trop_X_SS <- mu.link(Trop = 1,        CT = 0,       DD = 0,              IC = 0,       MI = 0,      PC = 0,        PP = 0,        SS = 1,       SP = 0, 
                        Trop_X_CT = 0,Trop_X_DD = 0,Trop_X_IC = 0, Trop_X_MI = 0,Trop_X_PC = 0,Trop_X_PP = 0, Trop_X_SS = 1, Trop_X_SP = 0
)


mu.Tem_X_SS <- mu.link(Trop = 0,        CT = 0,       DD = 0,              IC = 0,       MI = 0,      PC = 1,        PP = 0,        SS = 1,       SP = 0, 
                       Trop_X_CT = 0,Trop_X_DD = 0,Trop_X_IC = 0,Trop_X_MI = 0,Trop_X_PC = 0,Trop_X_PP = 0, Trop_X_SS = 0, Trop_X_SP = 0
)

dens( mu.Trop_X_SS, show.HPDI = .95, col="red",xlab = "N papers", main = "Sexual selection", xlim =c(1,25))
dens( mu.Tem_X_SS, show.HPDI = .95, add = T ) 



### Speciation 
mu.Trop_X_SP <- mu.link(Trop = 1,        CT = 0,       DD = 0,              IC = 0,       MI = 0,      PC = 0,        PP = 0,        SS = 0,       SP = 1, 
                        Trop_X_CT = 0,Trop_X_DD = 0,Trop_X_IC = 0, Trop_X_MI = 0,Trop_X_PC = 0,Trop_X_PP = 0, Trop_X_SS = 0, Trop_X_SP = 1
)


mu.Tem_X_SP <- mu.link(Trop = 0,        CT = 0,       DD = 0,              IC = 0,       MI = 0,      PC = 1,        PP = 0,        SS = 0,       SP = 1, 
                       Trop_X_CT = 0,Trop_X_DD = 0,Trop_X_IC = 0,Trop_X_MI = 0,Trop_X_PC = 0,Trop_X_PP = 0, Trop_X_SS = 0, Trop_X_SP = 0
)

dens( mu.Trop_X_SP, show.HPDI = .95, col="red",xlab = "N papers", main = "Speciation", xlim =c(1,20))
dens( mu.Tem_X_SP, show.HPDI = .95, add = T ) 





