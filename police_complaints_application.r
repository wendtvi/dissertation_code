#install.packages("staggered")
library(staggered)

setwd("C:/Users/vitor/Downloads")
load(file='pj_officer_level_balanced.rda')
dados=as.data.frame(pj_officer_level_balanced)
dados=dados[order(dados$uid,decreasing = FALSE),]

minha_base=data.frame(unique(dados$uid),0,0,0,0)
names(minha_base)=c("ID","COMPLAINTS_BEFORE","COMPLAINTS_AFTER","INTERVENTION_PERIOD","TREATED_GROUP")

c=1
k=0
for (i in unique(dados$uid)){
  k=k+1
  while (dados$uid[c]==i) {
    minha_base$INTERVENTION_PERIOD[k]=dados$first_trained[c]
    if(dados$period[c]<=dados$first_trained[c]) minha_base$COMPLAINTS_BEFORE[k]=minha_base$COMPLAINTS_BEFORE[k]+dados$complaints[c]
    if(dados$period[c]>dados$first_trained[c]) minha_base$COMPLAINTS_AFTER[k]=minha_base$COMPLAINTS_AFTER[k]+dados$complaints[c]
    c=c+1
  }
}

#escolho o período 28 como unico periodo de tratamento (quem recebeu tratamento em t = 28 é grupo tratado, quem recebeu após t=8 é grupo controle)
minha_base$TREATED_GROUP=as.numeric(minha_base$INTERVENTION_PERIOD==28)
minha_base=minha_base[minha_base$INTERVENTION_PERIOD>=28,]

#####################################
#############CIC ESTIMATOR###########
#####################################
Y_10=minha_base$COMPLAINTS_BEFORE[minha_base$TREATED_GROUP==1]
Y_11=minha_base$COMPLAINTS_AFTER[minha_base$TREATED_GROUP==1]
Y_00=minha_base$COMPLAINTS_BEFORE[minha_base$TREATED_GROUP==0]
Y_01=minha_base$COMPLAINTS_AFTER[minha_base$TREATED_GROUP==0]

#LB
F00_hat_LB=vector()
for (i in 1:length(Y_10)) { F00_hat_LB[i]=mean(Y_00<=Y_10[i]) }
F_inver_01_hat_LB=vector()
F01_hat_LB=vector()
for (i in 1:length(Y_01)) { F01_hat_LB[i]=mean(Y_01<=Y_01[i]) }
for (i in 1:length(F00_hat_LB)) { F_inver_01_hat_LB[i]=min(Y_01[F01_hat_LB>=F00_hat_LB[i]]) }
ATT_CIC_UB=mean(Y_11)-mean(F_inver_01_hat_LB)

#UB
F00_hat_UB=vector()
for (i in 1:length(Y_10)) { F00_hat_UB[i]=mean(Y_00<Y_10[i]) }
F_inver_01_hat_UB=vector()
F01_hat_UB=vector()
for (i in 1:length(Y_01)) { F01_hat_UB[i]=mean(Y_01<Y_01[i]) }
for (i in 1:length(F00_hat_UB)) { F_inver_01_hat_UB[i]=min(Y_01[F01_hat_UB>F00_hat_UB[i]]) }
ATT_CIC_LB=mean(Y_11)-mean(F_inver_01_hat_UB)



#####################################
#############MEU ESTIMATOR###########
#####################################
Y_10=minha_base$COMPLAINTS_BEFORE[minha_base$TREATED_GROUP==1]
Y_11=minha_base$COMPLAINTS_AFTER[minha_base$TREATED_GROUP==1]
Y_00=minha_base$COMPLAINTS_BEFORE[minha_base$TREATED_GROUP==0]
Y_01=minha_base$COMPLAINTS_AFTER[minha_base$TREATED_GROUP==0]

F00_hat_meu=ppois(Y_10, lambda = mean(Y_00))
F_inver_01_hat_meu=qpois(F00_hat_meu, lambda = mean(Y_01))

ATT_meu=mean(Y_11)-mean(F_inver_01_hat_meu)



#####################################
########WOOLDRIDGE ESTIMATOR#########
#####################################
Y_10=minha_base$COMPLAINTS_BEFORE[minha_base$TREATED_GROUP==1]
Y_11=minha_base$COMPLAINTS_AFTER[minha_base$TREATED_GROUP==1]
Y_00=minha_base$COMPLAINTS_BEFORE[minha_base$TREATED_GROUP==0]
Y_01=minha_base$COMPLAINTS_AFTER[minha_base$TREATED_GROUP==0]

ATT_WOOLD=mean(Y_11)-exp(log(mean(Y_10))+(log(mean(Y_01))-log(mean(Y_00))))


#################################################################
#resumo resultados
#CIC
ATT_CIC_UB
ATT_CIC_LB
(ATT_CIC_UB+ATT_CIC_LB)/2

#MEU
ATT_meu

#WOOLDRIDGE
ATT_WOOLD

