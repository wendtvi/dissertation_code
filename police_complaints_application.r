#install.packages("staggered")
library(staggered)

setwd("C:/Users/vitor/Downloads")
load(file='pj_officer_level_balanced.rda')
dados=as.data.frame(pj_officer_level_balanced)
dados=dados[dados$first_trained>13,]
dados=dados[dados$first_trained<72,]

matriz_estimadores_dinamicas=matrix(NA,length(unique(dados$first_trained))-1,14)
matriz_estimadores_dinamicas_media=matrix(NA,length(unique(dados$first_trained))-1,4)

l_matriz=0
max_period=max(dados$period)
for (p in sort(unique(dados$first_trained))[-1]){
  #  for (p in sort(unique(dados$first_trained))[2]){
  l_matriz=l_matriz+1
  t_interesse_1=sort(unique(dados$first_trained))[which(sort(unique(dados$first_trained))==p)-1]
  
  
  ##################################################
  ############EXPERIMENTO 2X2#######################
  ##################################################
  t_interesse=p
  
  c=1
  k=0
  i=0
  dados_temp=dados[dados$first_trained>t_interesse_1,]
  minha_base=data.frame(unique(dados_temp$uid),0,0,0,0)
  names(minha_base)=c("ID","complaints_BEFORE","complaints_AFTER","INTERVENTION_PERIOD","TREATED_GROUP")
  for (i in unique(dados_temp$uid)){
    k=k+1
    while (dados_temp$uid[c]==i && c<=length(dados_temp$uid)) {
      minha_base$INTERVENTION_PERIOD[k]=dados_temp$first_trained[c]
      if(dados_temp$period[c]<=t_interesse_1 && dados_temp$period[c]>t_interesse_1-5) minha_base$complaints_BEFORE[k]= minha_base$complaints_BEFORE[k]+dados_temp$complaints[c]
      if(dados_temp$period[c]>t_interesse && dados_temp$period[c]<=t_interesse+5) minha_base$complaints_AFTER[k]= minha_base$complaints_AFTER[k]+dados_temp$complaints[c]
      c=c+1
    }
  }
  
  #escolho o período 28 como unico periodo de tratamento (quem recebeu tratamento em t = 28 é grupo tratado, quem recebeu após t=28 é grupo controle)
  minha_base$TREATED_GROUP=as.numeric(minha_base$INTERVENTION_PERIOD==t_interesse)
  #minha_base=minha_base[minha_base$INTERVENTION_PERIOD>=t_interesse,]
  
  #####################################
  #############CIC ESTIMATOR###########
  #####################################
  Y_10=minha_base$complaints_BEFORE[minha_base$TREATED_GROUP==1]
  Y_11=minha_base$complaints_AFTER[minha_base$TREATED_GROUP==1]
  Y_00=minha_base$complaints_BEFORE[minha_base$TREATED_GROUP==0]
  Y_01=minha_base$complaints_AFTER[minha_base$TREATED_GROUP==0]
  
  #LB
  F00_hat_LB=vector()
  for (i in 1:length(Y_10)) { F00_hat_LB[i]=mean(Y_00<=Y_10[i]) }
  F_inver_01_hat_LB=vector()
  F01_hat_LB=vector()
  for (i in 1:length(Y_01)) { F01_hat_LB[i]=mean(Y_01<=Y_01[i]) }
  for (i in 1:length(F00_hat_LB)) { F_inver_01_hat_LB[i]=min(Y_01[F01_hat_LB>=F00_hat_LB[i]]) }
  ATT_CIC_LB=mean(Y_11)-mean(F_inver_01_hat_LB)*length(Y_10)/length(F_inver_01_hat_LB)
  
  #UB
  F00_hat_UB=vector()
  for (i in 1:length(Y_10)) { F00_hat_UB[i]=mean(Y_00<Y_10[i]) }
  F_inver_01_hat_UB=vector()
  F01_hat_UB=vector()
  for (i in 1:length(Y_01)) { F01_hat_UB[i]=mean(Y_01<Y_01[i]) }
  for (i in 1:length(F00_hat_UB)) { F_inver_01_hat_UB[i]=min(Y_01[F01_hat_UB>F00_hat_UB[i]]) }
  ATT_CIC_UB=mean(Y_11)-mean(F_inver_01_hat_UB)*length(Y_10)/length(F_inver_01_hat_UB)
  
  
  
  #####################################
  #############MEU ESTIMATOR###########
  #####################################
  Y_10=minha_base$complaints_BEFORE[minha_base$TREATED_GROUP==1]
  Y_11=minha_base$complaints_AFTER[minha_base$TREATED_GROUP==1]
  Y_00=minha_base$complaints_BEFORE[minha_base$TREATED_GROUP==0]
  Y_01=minha_base$complaints_AFTER[minha_base$TREATED_GROUP==0]
  
  pi_00=(var(Y_00)-sum(Y_00==0))/(var(Y_00)-sum(Y_00==0)^2-sum(Y_00==0))
  pi_01=(var(Y_01)-sum(Y_01==0))/(var(Y_01)-sum(Y_01==0)^2-sum(Y_01==0))
  pi_11=(var(Y_11)-sum(Y_11==0))/(var(Y_11)-sum(Y_11==0)^2-sum(Y_11==0))
  
  F00_hat_meu=ppois((Y_10), lambda = (1-pi_00)*mean(Y_00))
  F_inver_01_hat_meu=qpois((F00_hat_meu), lambda = (1-pi_01)*mean(Y_01))

  pi_F_inver=(var(F_inver_01_hat_meu)-sum(F_inver_01_hat_meu==0))/(var(F_inver_01_hat_meu)-sum(F_inver_01_hat_meu==0)^2-sum(F_inver_01_hat_meu==0))
  
  ATT_meu=(1-pi_11)*mean(Y_11)-(1-pi_F_inver)*mean(F_inver_01_hat_meu)
  
  
  
  #####################################
  ########WOOLDRIDGE ESTIMATOR#########
  #####################################
  Y_10=minha_base$complaints_BEFORE[minha_base$TREATED_GROUP==1]
  Y_11=minha_base$complaints_AFTER[minha_base$TREATED_GROUP==1]
  Y_00=minha_base$complaints_BEFORE[minha_base$TREATED_GROUP==0]
  Y_01=minha_base$complaints_AFTER[minha_base$TREATED_GROUP==0]
  
  ATT_WOOLD=mean(Y_11)-exp(log(mean(Y_10))+(log(mean(Y_01))-log(mean(Y_00))))
  
  
  #################################################################
  #resumo resultados
  #CIC
  matriz_estimadores_dinamicas[l_matriz,1]=ATT_CIC_UB
  matriz_estimadores_dinamicas[l_matriz,2]=ATT_CIC_LB
  
  #MEU
  matriz_estimadores_dinamicas[l_matriz,3]=ATT_meu
  
  #WOOLDRIDGE
  matriz_estimadores_dinamicas[l_matriz,4]=ATT_WOOLD
  
  matriz_estimadores_dinamicas[l_matriz,5]=length(Y_00)
  matriz_estimadores_dinamicas[l_matriz,6]=length(Y_10)
  matriz_estimadores_dinamicas[l_matriz,7]=mean(Y_00)
  matriz_estimadores_dinamicas[l_matriz,8]=mean(Y_01)
  matriz_estimadores_dinamicas[l_matriz,9]=mean(Y_10)
  matriz_estimadores_dinamicas[l_matriz,10]=mean(Y_11)
  matriz_estimadores_dinamicas[l_matriz,11]=mean(F_inver_01_hat_meu)
  matriz_estimadores_dinamicas[l_matriz,12]=mean(exp(log(mean(Y_10))+(log(mean(Y_01))-log(mean(Y_00)))))
  matriz_estimadores_dinamicas[l_matriz,13]=mean((F_inver_01_hat_LB))
  
}
matriz_estimadores_dinamicas[,14]=matriz_estimadores_dinamicas[,5]/(matriz_estimadores_dinamicas[,6]+matriz_estimadores_dinamicas[,5])
matriz_estimadores_dinamicas_temp=matriz_estimadores_dinamicas
matriz_estimadores_dinamicas=matriz_estimadores_dinamicas_temp[1:36,]
matriz_estimadores_dinamicas[matriz_estimadores_dinamicas[,1]==-Inf,1]=-10
matriz_estimadores_dinamicas[matriz_estimadores_dinamicas[,3]==-Inf,1]=-10
par(mfrow=c(1,1))
plot(y=matriz_estimadores_dinamicas[,3],x=sort(unique(dados$first_trained))[2:37], ylim = c(-2,1),pch=16,xlab = "Tempo",ylab = "ATT", main = "Janela temporal igual a 20")
points(y=matriz_estimadores_dinamicas[,4],x=sort(unique(dados$first_trained))[2:37],col="red",pch=17)
points(y=matriz_estimadores_dinamicas[,1],x=sort(unique(dados$first_trained))[2:37],col="blue",pch=24)
points(y=matriz_estimadores_dinamicas[,2],x=sort(unique(dados$first_trained))[2:37],col="blue",pch=25)

legend("topleft", legend=c("CICp", "Wooldridge", "CIC LI","CIC LS"),
       col=c("black","red", "blue","blue"),pch = c(16,17,24,25))















##################################################
############EXPERIMENTO STAGGERED#################
##################################################
#CS
mean_pre_treat_complaints=mean(dados$complaints[dados$period<dados$first_trained])
mean_pre_treat_complaints=mean(dados$complaints[dados$period<dados$first_trained])
mean_pre_treat_sustained=mean(dados$sustained[dados$period<dados$first_trained])

minha_base_staggered=dados
exposure_len_vector=72-minha_base_staggered$first_trained
hist(x = exposure_len_vector, main = "Comprimento do tempo de exposição",breaks = seq(0,60,5))

#PARALLEL TREND


staggered_cs(
  minha_base_staggered,
  i = "uid",
  t = "period",
  g = "first_trained",
  y = "force",
  estimand = "eventstudy",
  A_theta_list = NULL,
  A_0_list = NULL,
  eventTime = 0,
  return_full_vcv = FALSE,
  return_matrix_list = FALSE,
  compute_fisher = FALSE,
  num_fisher_permutations = 500,
  skip_data_check = FALSE
)/mean_pre_treat_force


#####################################
#############MEU ESTIMATOR###########
#####################################
#poisson
#complaints
vetor_ATT_treat_period=vector()
c=0
media_complaintsY11=vector()
media_complaintsY10=vector()
media_complaintsY01=vector()
media_complaintsY00=vector()
variancia_complaintsY11=vector()
variancia_complaintsY10=vector()
variancia_complaintsY01=vector()
variancia_complaintsY00=vector()
minha_base_staggered_temp=as.data.frame(minha_base_staggered)
for (i in sort(unique(minha_base_staggered$first_trained),decreasing = F)[-c(1:5,43:46)]){
  c=c+1
  Y_10=minha_base_staggered$complaints[(minha_base_staggered$period<i & minha_base_staggered$first_trained==i)]
  Y_11=minha_base_staggered$complaints[(minha_base_staggered$period>=i & minha_base_staggered$first_trained==i)]
  Y_00=minha_base_staggered$complaints[(minha_base_staggered$period<i & minha_base_staggered$first_trained>i)]
  Y_01=minha_base_staggered$complaints[(minha_base_staggered$period>=i & minha_base_staggered$first_trained>i & minha_base_staggered$period<minha_base_staggered$first_trained)]
  
  F00_hat_meu=pnbinom(Y_10, mu = mean(Y_00),size=length(Y_00))
  F_inver_01_hat_meu=qnbinom(F00_hat_meu, mu = mean(Y_01),size=length(Y_01))
  
  vetor_ATT_treat_period[c]=mean(Y_11)-mean(F_inver_01_hat_meu)
  
  media_complaintsY11[c]=mean(Y_11)
  media_complaintsY10[c]=mean(Y_10)
  media_complaintsY01[c]=mean(Y_01)
  media_complaintsY00[c]=mean(Y_00)
  variancia_complaintsY11[c]=var(Y_11)
  variancia_complaintsY10[c]=var(Y_10)
  variancia_complaintsY01[c]=var(Y_01)
  variancia_complaintsY00[c]=var(Y_00)
}

vetor_ATT_treat_period
#mean(vetor_ATT_treat_period[-((length(vetor_ATT_treat_period)-2):length(vetor_ATT_treat_period))])/mean_pre_treat_complaints
mean(vetor_ATT_treat_period,na.rm = T)/mean_pre_treat_complaints
(vetor_ATT_treat_period)/mean_pre_treat_complaints
plot(y=((vetor_ATT_treat_period)/mean_pre_treat_complaints),x=sort(unique(minha_base_staggered$first_trained),decreasing = F)[-c(1:5,43:46)],type='l')



#force
media_forceY11=vector()
media_forceY10=vector()
media_forceY01=vector()
media_forceY00=vector()
variancia_forceY11=vector()
variancia_forceY10=vector()
variancia_forceY01=vector()
variancia_forceY00=vector()
vetor_ATT_treat_period=vector()
c=0
for (i in sort(unique(minha_base_staggered$first_trained),decreasing = F)[-c(1:5,43:46)]){
  c=c+1
  Y_10=minha_base_staggered$force[(minha_base_staggered$period<i & minha_base_staggered$first_trained==i)]
  Y_11=minha_base_staggered$force[(minha_base_staggered$period>=i & minha_base_staggered$first_trained==i)]
  Y_00=minha_base_staggered$force[(minha_base_staggered$period<i & minha_base_staggered$first_trained>i)]
  Y_01=minha_base_staggered$force[(minha_base_staggered$period>=i & minha_base_staggered$first_trained>i & minha_base_staggered$period<minha_base_staggered$first_trained)]
  
  F00_hat_meu=pnbinom(Y_10, mu = mean(Y_00),size=length(Y_00))
  F_inver_01_hat_meu=qnbinom(F00_hat_meu, mu = mean(Y_01),size=length(Y_01))
  
  vetor_ATT_treat_period[c]=mean(Y_11)-mean(F_inver_01_hat_meu)
  
  media_forceY11[c]=mean(Y_11)
  media_forceY10[c]=mean(Y_10)
  media_forceY01[c]=mean(Y_01)
  media_forceY00[c]=mean(Y_00)
  variancia_forceY11[c]=var(Y_11)
  variancia_forceY10[c]=var(Y_10)
  variancia_forceY01[c]=var(Y_01)
  variancia_forceY00[c]=var(Y_00)
}

vetor_ATT_treat_period
#mean(vetor_ATT_treat_period[-((length(vetor_ATT_treat_period)-2):length(vetor_ATT_treat_period))])/mean_pre_treat_force
mean(vetor_ATT_treat_period,na.rm = T)/mean_pre_treat_force
(vetor_ATT_treat_period)/mean_pre_treat_force
plot(y=((vetor_ATT_treat_period)/mean_pre_treat_force),x=sort(unique(minha_base_staggered$first_trained),decreasing = F)[-c(1:5,43:46)],type='l')




#sustained
media_sustainedY11=vector()
media_sustainedY10=vector()
media_sustainedY01=vector()
media_sustainedY00=vector()
variancia_sustainedY11=vector()
variancia_sustainedY10=vector()
variancia_sustainedY01=vector()
variancia_sustainedY00=vector()
vetor_ATT_treat_period=vector()
c=0
for (i in sort(unique(minha_base_staggered$first_trained),decreasing = F)[-c(1:5,43:46)]){
  c=c+1
  Y_10=minha_base_staggered$sustained[(minha_base_staggered$period<i & minha_base_staggered$first_trained==i)]
  Y_11=minha_base_staggered$sustained[(minha_base_staggered$period>=i & minha_base_staggered$first_trained==i)]
  Y_00=minha_base_staggered$sustained[(minha_base_staggered$period<i & minha_base_staggered$first_trained>i)]
  Y_01=minha_base_staggered$sustained[(minha_base_staggered$period>=i & minha_base_staggered$first_trained>i & minha_base_staggered$period<minha_base_staggered$first_trained)]
  
  F00_hat_meu=pnbinom(Y_10, mu = mean(Y_00),size=length(Y_00))
  F_inver_01_hat_meu=qnbinom(F00_hat_meu, mu = mean(Y_01),size=length(Y_01))
  
  vetor_ATT_treat_period[c]=mean(Y_11)-mean(F_inver_01_hat_meu)
  
  media_sustainedY11[c]=mean(Y_11)
  media_sustainedY10[c]=mean(Y_10)
  media_sustainedY01[c]=mean(Y_01)
  media_sustainedY00[c]=mean(Y_00)
  variancia_sustainedY11[c]=var(Y_11)
  variancia_sustainedY10[c]=var(Y_10)
  variancia_sustainedY01[c]=var(Y_01)
  variancia_sustainedY00[c]=var(Y_00)
}

vetor_ATT_treat_period
#mean(vetor_ATT_treat_period[-((length(vetor_ATT_treat_period)-2):length(vetor_ATT_treat_period))])/mean_pre_treat_sustained
mean(vetor_ATT_treat_period,na.rm = T)/mean_pre_treat_sustained
(vetor_ATT_treat_period)/mean_pre_treat_sustained
plot(y=((vetor_ATT_treat_period)/mean_pre_treat_sustained),x=sort(unique(minha_base_staggered$first_trained),decreasing = F)[-c(1:5,43:46)],type='l')





#####################################
#############CIC ESTIMATOR###########
#####################################
c=0
vetor_ATT_treat_period_LB=vector()
vetor_ATT_treat_period_UB=vector()
for (i in sort(unique(minha_base_staggered$first_trained),decreasing = F)[-c(1:5,43:46)]){
  c=c+1
  Y_10=minha_base_staggered$complaints[(minha_base_staggered$period<i & minha_base_staggered$first_trained==i)]
  Y_11=minha_base_staggered$complaints[(minha_base_staggered$period>=i & minha_base_staggered$first_trained==i)]
  Y_00=minha_base_staggered$complaints[(minha_base_staggered$period<i & minha_base_staggered$first_trained>i)]
  Y_01=minha_base_staggered$complaints[(minha_base_staggered$period>=i & minha_base_staggered$first_trained>i & minha_base_staggered$period<minha_base_staggered$first_trained)]
  
  #LB
  F00_hat_LB=vector()
  for (j in 1:length(Y_10)) { F00_hat_LB[j]=mean(Y_00<=Y_10[j]) }
  F_inver_01_hat_LB=vector()
  F01_hat_LB=vector()
  for (j in 1:length(Y_01)) { F01_hat_LB[j]=mean(Y_01<=Y_01[j]) }
  for (j in 1:length(F00_hat_LB)) { F_inver_01_hat_LB[j]=min(Y_01[F01_hat_LB>=F00_hat_LB[j]]) }
  vetor_ATT_treat_period_LB[c]=mean(Y_11)-mean(F_inver_01_hat_LB)
  
  #UB
  F00_hat_UB=vector()
  for (j in 1:length(Y_10)) { F00_hat_UB[j]=mean(Y_00<Y_10[j]) }
  F_inver_01_hat_UB=vector()
  F01_hat_UB=vector()
  for (j in 1:length(Y_01)) { F01_hat_UB[j]=mean(Y_01<Y_01[j]) }
  for (j in 1:length(F00_hat_UB)) { F_inver_01_hat_UB[j]=min(Y_01[F01_hat_UB>F00_hat_UB[j]]) }
  vetor_ATT_treat_period_UB[c]=mean(Y_11)-mean(F_inver_01_hat_UB)
  
}
vetor_ATT_treat_period_UB_temp=vector()
vetor_ATT_treat_period_UB_temp=vetor_ATT_treat_period_UB
vetor_ATT_treat_period_UB_temp[vetor_ATT_treat_period_UB_temp < (-10)]=-0.9761905
vetor_ATT_treat_period_LB
#mean(vetor_ATT_treat_period[-((length(vetor_ATT_treat_period)-2):length(vetor_ATT_treat_period))])/mean_pre_treat_sustained
mean(vetor_ATT_treat_period_LB,na.rm = T)/mean_pre_treat_complaints
(vetor_ATT_treat_period_LB)/mean_pre_treat_complaints
plot(y=((vetor_ATT_treat_period_LB)),x=sort(unique(minha_base_staggered$first_trained),decreasing = F)[-c(1:5,43:46)],type='l',ylim = c(-1.1,0),ylab = "",xlab = "")
lines(y=((vetor_ATT_treat_period_UB_temp)),x=sort(unique(minha_base_staggered$first_trained),decreasing = F)[-c(1:5,43:46)],col="red")

par(mfrow=c(1,2))
plot(y=((vetor_ATT_treat_period_LB)),x=sort(unique(minha_base_staggered$first_trained),decreasing = F)[-c(1:5,43:46)],type='l',ylab = "",xlab = "")
plot(y=((vetor_ATT_treat_period_UB_temp)),x=sort(unique(minha_base_staggered$first_trained),decreasing = F)[-c(1:5,43:46)],type='l',col="red",ylab = "",xlab = "")



########################################
##############TESTES####################
########################################
#Teste independencia (idade do policial)
plot(y=minha_base_staggered$birth_year, x=minha_base_staggered$first_trained)
vetor_medias_ano=vector()
for (k in sort(unique(minha_base_staggered$first_trained),decreasing = F)){
  vetor_medias_ano[k]=mean(minha_base_staggered$birth_year[minha_base_staggered$first_trained==k])
} 
points(vetor_medias_ano,col="red")


#Teste distribuição poisson
par(mfrow=c(2,2))
plot((media_complaintsY11/variancia_complaintsY11),type='l',ylab = "")
plot((media_complaintsY10/variancia_complaintsY10),type='l',ylab = "")
plot((media_complaintsY01/variancia_complaintsY01),type='l',ylab = "")
plot((media_complaintsY00/variancia_complaintsY00),type='l',ylab = "")



#install.packages("vcd")
library(vcd) ## loading vcd package

i=31
Y_10=minha_base_staggered$complaints[(minha_base_staggered$period<i & minha_base_staggered$first_trained==i)]
Y_11=minha_base_staggered$complaints[(minha_base_staggered$period>=i & minha_base_staggered$first_trained==i)]
Y_00=minha_base_staggered$complaints[(minha_base_staggered$period<i & minha_base_staggered$first_trained>i)]
Y_01=minha_base_staggered$complaints[(minha_base_staggered$period>=i & minha_base_staggered$first_trained>i & minha_base_staggered$period<minha_base_staggered$first_trained)]

y=Y_11 
gf = goodfit(y,type= "nbinomial",method= "ML")
plot(gf,main="Count data vs Poisson distribution")
summary(gf)
y=Y_10 
gf = goodfit(y,type= "nbinomial",method= "ML")
plot(gf,main="Count data vs Poisson distribution")
summary(gf)
y=Y_01 
gf = goodfit(y,type= "nbinomial",method= "ML")
plot(gf,main="Count data vs Poisson distribution")
summary(gf)
y=Y_00 
gf = goodfit(y,type= "nbinomial",method= "ML")
plot(gf,main="Count data vs Poisson distribution")
summary(gf)


#Teste invariancia no tempo
i=22

Y_10=minha_base_staggered$complaints[(minha_base_staggered$period<i & minha_base_staggered$first_trained==i)]
Y_11=minha_base_staggered$complaints[(minha_base_staggered$period>=i & minha_base_staggered$first_trained==i)]
Y_00=minha_base_staggered$complaints[(minha_base_staggered$period<i & minha_base_staggered$first_trained>i)]
Y_01=minha_base_staggered$complaints[(minha_base_staggered$period>=i & minha_base_staggered$first_trained>i & minha_base_staggered$period<minha_base_staggered$first_trained)]
F00_hat_meu=pnbinom(Y_10, mu = mean(Y_00),size=length(Y_00))
F_inver_01_hat_meu=qnbinom(F00_hat_meu, mu = mean(Y_01),size=length(Y_01))

par(mfrow=c(2,2))
y=F_inver_01_hat_meu 
gf = goodfit(y,type= "nbinomial",method= "ML")
plot(gf,main="Count data vs NBinomial distribution")
summary(gf)
y=Y_10 
gf = goodfit(y,type= "nbinomial",method= "ML")
plot(gf,main="Count data vs NBinomial distribution")
summary(gf)
y=Y_01
gf = goodfit(y,type= "nbinomial",method= "ML")
plot(gf,main="Count data vs NBinomial distribution")
summary(gf)
y=Y_00
gf = goodfit(y,type= "nbinomial",method= "ML")
plot(gf,main="Count data vs NBinomial distribution")
summary(gf)
