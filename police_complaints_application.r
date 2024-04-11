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
  
  pi_00=(var(Y_00)-mean(Y_00))/(var(Y_00)+mean(Y_00)^2-mean(Y_00))
  pi_01=(var(Y_01)-mean(Y_01))/(var(Y_01)+mean(Y_01)^2-mean(Y_01))
  pi_11=(var(Y_11)-mean(Y_11))/(var(Y_11)+mean(Y_11)^2-mean(Y_11))
  
  F00_hat_meu=ppois(sum(Y_10), lambda = (1-pi_00)*sum(Y_00))
  F_inver_01_hat_meu=qpois((F00_hat_meu), lambda = (1-pi_01)*sum(Y_01))

  pi_F_inver=(var(F_inver_01_hat_meu)-mean(F_inver_01_hat_meu))/(var(F_inver_01_hat_meu)+mean(F_inver_01_hat_meu)^2-mean(F_inver_01_hat_meu))
    
  ATT_meu=(1-pi_11)*mean(Y_11)-(1-pi_F_inver)*sum(F_inver_01_hat_meu)/length(F_inver_01_hat_meu)
  
  
  
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

