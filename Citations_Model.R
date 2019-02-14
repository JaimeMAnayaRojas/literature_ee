####  
library(rethinking)
library(data.table)
rm(list=ls(all=TRUE))


data<-fread("./data/AllData.csv", header=T, colClasses="character", na.strings = "")
names(data)

data$TOPIC = factor(data$TOPIC)
levels(data$TOPIC)


levels(data$TOPIC) <- c("AD", "CT", "DD", "IC", "MI", "PC", "PP", "SS", "SP")

data <- data[-which(data$TOPIC == "PP"), ]

data$TAXON = factor(data$TAXON)
levels(data$TAXON)

data$TOPIC = factor(data$TOPIC)

data$ID <- 1:dim(data)[1]


#get the number of publications per topic, region, taxon
mod.mat <- as.data.frame(model.matrix(~ REGION*TOPIC*TAXON + PERIOD + ID, data))
data2 <- data[mod.mat$ID, ]
mod.mat <- as.data.frame(model.matrix(~ REGION*TOPIC + PERIOD, data2))

names(mod.mat) <- c("(Intercept)","TROP","CT","DD","IC","MI", "PC","SS","SP","P2","P3",
                    "P4","P5","P6","P7","TROP_X_CT","TROP_X_DD", "TROP_X_IC","TROP_X_MI","TROP_X_PC","TROP_X_SS","TROP_X_SP")




d <-as.data.frame(cbind(data2[,c("CITATIONS","TAXON") ], mod.mat[,-1]))
levels(d$TAXON)<- 1: length(levels(d$TAXON))
d$TAXON <- as.numeric(d$TAXON)
names(d)[1] <- "Cite"
str(d)
d$Cite <- as.numeric(d$Cite)

names(data2)
# get the bayesian model 
glimmer(CITATIONS ~ REGION*TOPIC + PERIOD + (1|TAXON) ,  data2)


####### Full model
model.1 <- map2stan(
  
  alist(
    Cite ~ dnorm( mu , sigma ),
    mu <- Intercept + b_TROP*TROP +
      b_CT*CT + b_DD*DD + b_IC*IC + b_MI*MI +
      b_PC*PC +  b_SS*SS + b_SP*SP +
      b_P2*P2 + b_P3*P3 + b_P4*P4 + b_P5*P5 +
      b_P6*P6 + b_P7*P7 + b_TROP_X_CT*TROP_X_CT +
      b_TROP_X_DD*TROP_X_DD +
      b_TROP_X_IC*TROP_X_IC +
      b_TROP_X_MI*TROP_X_MI +
      b_TROP_X_PC*TROP_X_PC +
      b_TROP_X_SS*TROP_X_SS +
      b_TROP_X_SP*TROP_X_SP +
      v_Intercept[TAXON],
    Intercept ~ dnorm(0,10),
    b_TROP ~ dnorm(0,10),
    b_CT ~ dnorm(0,10),
    b_DD ~ dnorm(0,10),
    b_IC ~ dnorm(0,10),
    b_MI ~ dnorm(0,10),
    b_PC ~ dnorm(0,10),
    b_SS ~ dnorm(0,10),
    b_SP ~ dnorm(0,10),
    b_P2 ~ dnorm(0,10),
    b_P3 ~ dnorm(0,10),
    b_P4 ~ dnorm(0,10),
    b_P5 ~ dnorm(0,10),
    b_P6 ~ dnorm(0,10),
    b_P7 ~ dnorm(0,10),
    b_TROP_X_CT ~ dnorm(0,10),
    b_TROP_X_DD ~ dnorm(0,10),
    b_TROP_X_IC ~ dnorm(0,10),
    b_TROP_X_MI ~ dnorm(0,10),
    b_TROP_X_PC ~ dnorm(0,10),
    b_TROP_X_SS ~ dnorm(0,10),
    b_TROP_X_SP ~ dnorm(0,10),
    v_Intercept[TAXON] ~ dnorm(0,sigma_TAXON),
    sigma_TAXON ~ dcauchy(0,2),
    sigma ~ dcauchy(0,2)
  ), data = d, chains = 4, iter = 5000, WAIC = T
  
  
)


precis(model.1, depth = 2, prob = .95, digits = 5)
## 






model.2 <- map2stan(
  
  alist(
    Cite ~ dnorm( mu , sigma ),
    mu <- Intercept + b_TROP*TROP +
      b_CT*CT + b_DD*DD + b_IC*IC + b_MI*MI +
      b_PC*PC + b_SS*SS + b_SP*SP +
      b_P2*P2 + b_P3*P3 + b_P4*P4 + b_P5*P5 +
      b_P6*P6 + b_P7*P7 +
      v_Intercept[TAXON],
    Intercept ~ dnorm(0,10),
    b_TROP ~ dnorm(0,10),
    b_CT ~ dnorm(0,10),
    b_DD ~ dnorm(0,10),
    b_IC ~ dnorm(0,10),
    b_MI ~ dnorm(0,10),
    b_PC ~ dnorm(0,10),
    b_SS ~ dnorm(0,10),
    b_SP ~ dnorm(0,10),
    b_P2 ~ dnorm(0,10),
    b_P3 ~ dnorm(0,10),
    b_P4 ~ dnorm(0,10),
    b_P5 ~ dnorm(0,10),
    b_P6 ~ dnorm(0,10),
    b_P7 ~ dnorm(0,10),
    v_Intercept[TAXON] ~ dnorm(0,sigma_TAXON),
    sigma_TAXON ~ dcauchy(0,2),
    sigma ~ dcauchy(0,2)
  ), data = d, chains = 4, iter = 5000, WAIC = T
  
  
)

#check if the model converges

# compare the models, is the interaction term important
compare(model.1, model.2)



# get the parameters
(sumStat <- precis(model.1, depth = 2, prob = .95, digits = 5))
sumStat[sumStat[,3]>0,]
sumStat[sumStat[,4]<0,]

write.csv(sumStat, "Parameters_CITATIONS.csv")

saveRDS(model.1, "data/Ciation_model_1.rds")
#model.1 <- readRDS("data/Ciation_model_1.rds")


# Get the posterior probabilities and plot the values for each treatment combination
pos <- extract.samples(model.1)

graphics.off()



mu.link <- function(TROP, CT, DD, IC, MI, PC, SS, SP, 
                    TROP_X_CT, TROP_X_DD, TROP_X_IC, TROP_X_MI, TROP_X_PC,  TROP_X_SS, TROP_X_SP
                    
){
  
  Y.est <- with(pos, Intercept + b_TROP*TROP + b_CT*CT + b_DD*DD + b_IC*IC + 
                  b_MI*MI + b_PC*PC +  + b_SS*SS+ b_SP*SP + 
                  b_TROP_X_CT*TROP_X_CT + b_TROP_X_DD*TROP_X_DD + b_TROP_X_IC*TROP_X_IC + 
                  b_TROP_X_MI*TROP_X_MI + b_TROP_X_PC*TROP_X_PC +  
                  b_TROP_X_SS*TROP_X_SS+ b_TROP_X_SP*TROP_X_SP 
  )
  return(exp(Y.est))
  
}




# Adaptation x region  
mu.TROP_X_AD <- mu.link(TROP = 1, CT = 0,DD = 0,IC = 0,MI = 0,PC = 0,SS = 0,SP = 0, 
                        TROP_X_CT = 0,TROP_X_DD = 0,TROP_X_IC = 0,TROP_X_MI = 0,TROP_X_PC = 0,
                        TROP_X_SS = 0, TROP_X_SP = 0
)

mu.Tem_X_AD <- mu.link(TROP = 0, CT = 0,DD = 0,IC = 0,MI = 0,PC = 0,SS = 0,SP = 0, 
                       TROP_X_CT = 0,TROP_X_DD = 0,TROP_X_IC = 0,TROP_X_MI = 0,TROP_X_PC = 0,
                       TROP_X_SS = 0, TROP_X_SP = 0
)




# Climate Tolerance x region  

mod.mat[which(mod.mat$TROP == 1 & mod.mat$CT ==1),]

mu.TROP_X_CT <- mu.link(TROP = 1,        CT = 1,       DD = 0,              IC = 0,       MI = 0,      PC = 0,          SS = 0,       SP = 0, 
                        TROP_X_CT = 1,TROP_X_DD = 0,TROP_X_IC = 0,TROP_X_MI = 0,TROP_X_PC = 0, TROP_X_SS = 0, TROP_X_SP = 0
)

mod.mat[which(mod.mat$TROP == 0 & mod.mat$CT ==1),]
mu.Tem_X_CT <- mu.link(TROP = 0,        CT = 1,       DD = 0,              IC = 0,       MI = 0,      PC = 0,      SS = 0,       SP = 0, 
                       TROP_X_CT = 0,TROP_X_DD = 0,TROP_X_IC = 0,TROP_X_MI = 0,TROP_X_PC = 0, TROP_X_SS = 0, TROP_X_SP = 0
)




# density x region  

mod.mat[which(mod.mat$TROP == 1 & mod.mat$DD ==1),]

mu.TROP_X_DD <- mu.link(TROP = 1,        CT = 0,       DD = 1,              IC = 0,       MI = 0,      PC = 0,      SS = 0,       SP = 0, 
                        TROP_X_CT = 1,TROP_X_DD = 1,TROP_X_IC = 0,TROP_X_MI = 0,TROP_X_PC = 0, TROP_X_SS = 0, TROP_X_SP = 0
)

mod.mat[which(mod.mat$TROP == 0 & mod.mat$DD ==1),]
mu.Tem_X_DD <- mu.link(TROP = 0,        CT = 0,       DD = 1,              IC = 0,       MI = 0,      PC = 0,       SS = 0,       SP = 0, 
                       TROP_X_CT = 0,TROP_X_DD = 0,TROP_X_IC = 0,TROP_X_MI = 0,TROP_X_PC = 0, TROP_X_SS = 0, TROP_X_SP = 0
)






# Interspecific comp x region  

mod.mat[which(mod.mat$TROP == 1 & mod.mat$IC ==1),]

mu.TROP_X_IC <- mu.link(TROP = 1,        CT = 0,       DD = 0,              IC = 1,       MI = 0,      PC = 0,     SS = 0,       SP = 0, 
                        TROP_X_CT = 1,TROP_X_DD = 1,TROP_X_IC = 1, TROP_X_MI = 0,TROP_X_PC = 0,TROP_X_SS = 0, TROP_X_SP = 0
)

mod.mat[which(mod.mat$TROP == 0 & mod.mat$IC ==1),]
mu.Tem_X_IC <- mu.link(TROP = 0,        CT = 0,       DD = 0,              IC = 1,       MI = 0,      PC = 0,        SS = 0,       SP = 0, 
                       TROP_X_CT = 0,TROP_X_DD = 0,TROP_X_IC = 0,TROP_X_MI = 0,TROP_X_PC = 0, TROP_X_SS = 0, TROP_X_SP = 0
)





# Mimicry x region  


mu.TROP_X_MI <- mu.link(TROP = 1,        CT = 0,       DD = 0,              IC = 0,       MI = 1,      PC = 0,        SS = 0,       SP = 0, 
                        TROP_X_CT = 0,TROP_X_DD = 0,TROP_X_IC = 0, TROP_X_MI = 1,TROP_X_PC = 0, TROP_X_SS = 0, TROP_X_SP = 0
)


mu.Tem_X_MI <- mu.link(TROP = 0,        CT = 0,       DD = 0,              IC = 0,       MI = 1,      PC = 0,       SS = 0,       SP = 0, 
                       TROP_X_CT = 0,TROP_X_DD = 0,TROP_X_IC = 0,TROP_X_MI = 0,TROP_X_PC = 0, TROP_X_SS = 0, TROP_X_SP = 0
)





### PARENTAL CARE
mu.TROP_X_PC <- mu.link(TROP = 1,        CT = 0,       DD = 0,              IC = 0,       MI = 0,      PC = 1,         SS = 0,       SP = 0, 
                        TROP_X_CT = 0,TROP_X_DD = 0,TROP_X_IC = 0, TROP_X_MI = 0,TROP_X_PC = 1,TROP_X_SS = 0, TROP_X_SP = 0
)


mu.Tem_X_PC <- mu.link(TROP = 0,        CT = 0,       DD = 0,              IC = 0,       MI = 0,      PC = 1,   SS = 0,   SP = 0, 
                       TROP_X_CT = 0,TROP_X_DD = 0,TROP_X_IC = 0,TROP_X_MI = 0,TROP_X_PC = 0, TROP_X_SS = 0, TROP_X_SP = 0
)



### Sexual selection 
mu.TROP_X_SS <- mu.link(TROP = 1,        CT = 0,       DD = 0,              IC = 0,       MI = 0,      PC = 0,    SS = 1,       SP = 0, 
                        TROP_X_CT = 0,TROP_X_DD = 0,TROP_X_IC = 0, TROP_X_MI = 0,TROP_X_PC = 0,TROP_X_SS = 1, TROP_X_SP = 0
)


mu.Tem_X_SS <- mu.link(TROP = 0,        CT = 0,       DD = 0,              IC = 0,       MI = 0, PC = 0, SS = 1,       SP = 0, 
                       TROP_X_CT = 0,TROP_X_DD = 0,TROP_X_IC = 0,TROP_X_MI = 0,TROP_X_PC = 0, TROP_X_SS = 0, TROP_X_SP = 0
)




### Speciation 
mu.TROP_X_SP <- mu.link(TROP = 1,        CT = 0,       DD = 0,              IC = 0,       MI = 0, PC = 0,  SS = 0, SP = 1, 
                        TROP_X_CT = 0,TROP_X_DD = 0,TROP_X_IC = 0, TROP_X_MI = 0,TROP_X_PC = 0, TROP_X_SS = 0, TROP_X_SP = 1
)


mu.Tem_X_SP <- mu.link(TROP = 0,        CT = 0,   DD = 0,  IC = 0, MI = 0,  PC = 0, SS = 0,       SP = 1, 
                       TROP_X_CT = 0,TROP_X_DD = 0,TROP_X_IC = 0,TROP_X_MI = 0,TROP_X_PC = 0,TROP_X_SS = 0, TROP_X_SP = 0
)







######## Plot the distribution of citations

svg(filename = "plots/Citations region by topic.svg", width = 697, height = 10.36 )
op<-par(mfrow=c(3,3),mar=c(4,2,2,2))

dens( mu.Tem_X_AD, show.HPDI = .95,  xlab = "Citation rate", main = "Adaptation", xlim=c(20,60))
dens( mu.TROP_X_AD, show.HPDI = .95,col="red",   xlim= c(0,60), add = T)
PerTembia <- paste(round(100 * length(which(mu.TROP_X_AD - mu.Tem_X_AD < 0))/ length(mu.Tem_X_AD),2), "%", sep = "")
legend("topright", PerTembia, bty = "n")



dens( mu.Tem_X_CT, show.HPDI = .95,  xlab = "Citation rate", main = "Climate tolerance", xlim = c(25, 75))
dens( mu.TROP_X_CT, show.HPDI = .95,col="red", add = T) 

PerTembia <- paste(round(100 * length(which(mu.TROP_X_CT - mu.Tem_X_CT < 0))/ length(mu.Tem_X_AD),2), "%", sep = "")
legend("topright", PerTembia, bty = "n")





dens( mu.Tem_X_DD, show.HPDI = .95, xlab = "Citation rate", main = "Density dependence", xlim= c(15, 75)) 
dens( mu.TROP_X_DD, show.HPDI = .95, col="red", add = T)
PerTembia <- paste(round(100 * length(which(mu.TROP_X_DD - mu.Tem_X_DD < 0))/ length(mu.Tem_X_DD),2), "%", sep = "")
legend("topright", PerTembia, bty = "n")



dens( mu.Tem_X_IC, show.HPDI = .95, xlab = "Citation rate", main = "Interspecific competition", xlim=c(10,70)) 
dens( mu.TROP_X_IC, show.HPDI = .95, col="red",add = T)
PerTembia <- paste(round(100 * length(which(mu.TROP_X_IC - mu.Tem_X_IC < 0))/ length(mu.Tem_X_AD),2), "%", sep = "")
legend("topright", PerTembia, bty = "n")



dens( mu.Tem_X_MI, show.HPDI = .95, xlab = "Citation rate", main = "Mimicry" )
dens( mu.TROP_X_MI, show.HPDI = .95, col="red",add = T)
PerTembia <- paste(round(100 * length(which(mu.TROP_X_MI - mu.Tem_X_MI < 0))/ length(mu.Tem_X_MI),2), "%", sep = "")
legend("topright", PerTembia, bty = "n")





dens( mu.Tem_X_PC, show.HPDI = .95, xlab = "Citation rate", main = "Parental care", xlim=c(5,50)) 
dens( mu.TROP_X_PC, show.HPDI = .95, col="red", add = T)
PerTembia <- paste(round(100 * length(which(mu.TROP_X_PC - mu.Tem_X_PC < 0))/ length(mu.Tem_X_PC),2), "%", sep = "")
legend("topright", PerTembia, bty = "n")


dens( mu.Tem_X_SS, show.HPDI = .95, xlab = "Citation rate", main = "Sexual selection")
dens(mu.TROP_X_SS, show.HPDI = .95, col="red" , add = T ) 
PerTembia <- paste(round(100 * length(which(mu.TROP_X_SS - mu.Tem_X_SS < 0))/ length(mu.Tem_X_SS),2), "%", sep = "")
legend("topright", PerTembia, bty = "n")


dens( mu.Tem_X_SP , show.HPDI = .95,xlab = "Citation rate", main = "Speciation", xlim=c(10,60), ylim=c(0,.09))
dens(mu.TROP_X_SP, col="red", show.HPDI = .95, add = T ) 

PerTembia <- paste(round(100 * length(which(mu.TROP_X_SP - mu.Tem_X_SP < 0))/ length(mu.Tem_X_SP),2), "%", sep = "")
legend("topright", PerTembia, bty = "n")

graphics.off()



