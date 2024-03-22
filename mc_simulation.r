##############################################
###############MCMC com heatmaps##############
##############################################
#Parâmetro de interesse heatmap = escala da distribuição do termo de erro do modelo latente
vetor_heatmap=vector()
vetor_parametro1=vector()
vetor_parametro2=vector()

cc=0
for (escala in seq(0,1,0.2)){
  for (forma in seq(1,2,0.2)){
    cc=cc+1
    #install.packages("devtools")
    #library(devtools)
    #install.packages("moments")
    #library(moments)
    
    #source_url("https://github.com/wendtvi/dissertation_code/blob/main/continous_nonlinear.r?raw=TRUE")
    N=100
    n1=500
    full_return=mc_function(N=N,parametro1=escala,parametro2=forma)
    WOOLmatriz_resultados=full_return[[1]]
    matriz_resultados=full_return[[2]]
    matriz_var_est_4=full_return[[3]]
    matriz_var_est_4_assint=matrix(NA,nrow = N,ncol = 2)
    matriz_var_est_4_assint=matriz_var_est_4[,1:2]^2*matriz_var_est_4[,3]/n1
    sqrt(mean(matriz_var_est_4[1,]))
    matriz_var_est_5=full_return[[4]]
    matriz_var_est_5_assint=matrix(NA,nrow = N,ncol = 2)
    matriz_var_est_5_assint=matriz_var_est_5[,1:2]^2*matriz_var_est_5[,3]/n1
    sqrt(mean(matriz_var_est_5[1,]))
    matriz_var_est_6=full_return[[5]]
    matriz_var_est_6_assint=matrix(NA,nrow = N,ncol = 2)
    matriz_var_est_6_assint=matriz_var_est_6[,1:2]^2*matriz_var_est_6[,3]/n1
    sqrt(mean(matriz_var_est_6[1,]))
    
    WOOLvetor_mc_resultados_vies=vector()
    WOOLvetor_mc_resultados_sd=vector()
    WOOLvetor_mc_resultados_media_pop=vector()
    WOOLvetor_mc_resultados_media_est=vector()
    WOOLvetor_mc_resultados_sd_vies=vector()
    WOOLvetor_mc_resultados_mcsd=vector()
    
    WOOLxvetor_mc_resultados_vies=vector()
    WOOLxvetor_mc_resultados_sd=vector()
    WOOLxvetor_mc_resultados_media_pop=vector()
    WOOLxvetor_mc_resultados_media_est=vector()
    WOOLxvetor_mc_resultados_sd_vies=vector()
    WOOLxvetor_mc_resultados_mcsd=vector()
    
    WOOLyvetor_mc_resultados_vies=vector()
    WOOLyvetor_mc_resultados_sd=vector()
    WOOLyvetor_mc_resultados_media_pop=vector()
    WOOLyvetor_mc_resultados_media_est=vector()
    WOOLyvetor_mc_resultados_sd_vies=vector()
    WOOLyvetor_mc_resultados_mcsd=vector()
    
    simulation_residuals=matrix(0,nrow = N,ncol = 2)
    
    c=0
    for (k in seq(1,ncol(WOOLmatriz_resultados)/3,2)){
      c=c+1
      WOOLvetor_mc_resultados_vies[c]=mean(WOOLmatriz_resultados[,k])-mean(WOOLmatriz_resultados[,k+1])
      WOOLvetor_mc_resultados_sd_vies[c]=mean(abs(WOOLmatriz_resultados[,k]-(WOOLmatriz_resultados[,k+1])))
      WOOLvetor_mc_resultados_mcsd[c]=var(WOOLmatriz_resultados[,k+1])*sqrt(((WOOLmatriz_resultados[,k+1])-1)/N)
      WOOLvetor_mc_resultados_sd[c]=sd(WOOLmatriz_resultados[,k+1])
      WOOLvetor_mc_resultados_media_pop[c]=mean(WOOLmatriz_resultados[,k])
      WOOLvetor_mc_resultados_media_est[c]=mean(WOOLmatriz_resultados[,k+1])
    }
    
    c=0
    for (k in seq(ncol(WOOLmatriz_resultados)/3+1,ncol(WOOLmatriz_resultados)*2/3,2)){
      c=c+1
      WOOLxvetor_mc_resultados_vies[c]=mean(WOOLmatriz_resultados[,k])-mean(WOOLmatriz_resultados[,k+1])
      WOOLxvetor_mc_resultados_sd_vies[c]=mean(abs(WOOLmatriz_resultados[,k]-(WOOLmatriz_resultados[,k+1])))
      WOOLxvetor_mc_resultados_mcsd[c]=var(WOOLmatriz_resultados[,k+1])*sqrt(((WOOLmatriz_resultados[,k+1])-1)/N)
      WOOLxvetor_mc_resultados_sd[c]=sd(WOOLmatriz_resultados[,k+1])
      WOOLxvetor_mc_resultados_media_pop[c]=mean(WOOLmatriz_resultados[,k])
      WOOLxvetor_mc_resultados_media_est[c]=mean(WOOLmatriz_resultados[,k+1])
    }
    
    
    c=0
    for (k in seq(ncol(WOOLmatriz_resultados)*2/3+1,ncol(WOOLmatriz_resultados),2)){
      c=c+1
      WOOLyvetor_mc_resultados_vies[c]=mean(WOOLmatriz_resultados[,k])-mean(WOOLmatriz_resultados[,k+1])
      WOOLyvetor_mc_resultados_sd_vies[c]=mean(abs(WOOLmatriz_resultados[,k]-(WOOLmatriz_resultados[,k+1])))
      WOOLyvetor_mc_resultados_mcsd[c]=var(WOOLmatriz_resultados[,k+1])*sqrt(((WOOLmatriz_resultados[,k+1])-1)/N)
      WOOLyvetor_mc_resultados_sd[c]=sd(WOOLmatriz_resultados[,k+1])
      WOOLyvetor_mc_resultados_media_pop[c]=mean(WOOLmatriz_resultados[,k])
      WOOLyvetor_mc_resultados_media_est[c]=mean(WOOLmatriz_resultados[,k+1])
    }
    
    
    vetor_mc_resultados_vies=vector()
    vetor_mc_resultados_sd=vector()
    vetor_mc_resultados_media_pop=vector()
    vetor_mc_resultados_media_est=vector()
    vetor_mc_resultados_sd_vies=vector()
    vetor_mc_resultados_mcsd=vector()
    
    xvetor_mc_resultados_vies=vector()
    xvetor_mc_resultados_sd=vector()
    xvetor_mc_resultados_media_pop=vector()
    xvetor_mc_resultados_media_est=vector()
    xvetor_mc_resultados_sd_vies=vector()
    xvetor_mc_resultados_mcsd=vector()
    
    yvetor_mc_resultados_vies=vector()
    yvetor_mc_resultados_sd=vector()
    yvetor_mc_resultados_media_pop=vector()
    yvetor_mc_resultados_media_est=vector()
    yvetor_mc_resultados_sd_vies=vector()
    yvetor_mc_resultados_mcsd=vector()
    
    c=0
    for (k in seq(1,ncol(matriz_resultados)/3,2)){
      c=c+1
      vetor_mc_resultados_vies[c]=mean(matriz_resultados[,k])-mean(matriz_resultados[,k+1])
      vetor_mc_resultados_sd_vies[c]=mean(abs(matriz_resultados[,k]-(matriz_resultados[,k+1])))
      vetor_mc_resultados_mcsd[c]=var(matriz_resultados[,k+1])*sqrt(((matriz_resultados[,k+1])-1)/N)
      vetor_mc_resultados_sd[c]=sd(matriz_resultados[,k+1])
      vetor_mc_resultados_media_pop[c]=mean(matriz_resultados[,k])
      vetor_mc_resultados_media_est[c]=mean(matriz_resultados[,k+1])
    }
    
    c=0
    for (k in seq(ncol(matriz_resultados)/3+1,ncol(matriz_resultados)*2/3,2)){
      c=c+1
      xvetor_mc_resultados_vies[c]=mean(matriz_resultados[,k])-mean(matriz_resultados[,k+1])
      xvetor_mc_resultados_sd_vies[c]=mean(abs(matriz_resultados[,k]-(matriz_resultados[,k+1])))
      xvetor_mc_resultados_mcsd[c]=var(matriz_resultados[,k+1])*sqrt(((matriz_resultados[,k+1])-1)/N)
      xvetor_mc_resultados_sd[c]=sd(matriz_resultados[,k+1])
      xvetor_mc_resultados_media_pop[c]=mean(matriz_resultados[,k])
      xvetor_mc_resultados_media_est[c]=mean(matriz_resultados[,k+1])
    }
    
    
    c=0
    for (k in seq(ncol(matriz_resultados)*2/3+1,ncol(matriz_resultados),2)){
      c=c+1
      yvetor_mc_resultados_vies[c]=mean(matriz_resultados[,k])-mean(matriz_resultados[,k+1])
      yvetor_mc_resultados_sd_vies[c]=mean(abs(matriz_resultados[,k]-(matriz_resultados[,k+1])))
      yvetor_mc_resultados_mcsd[c]=var(matriz_resultados[,k+1])*sqrt(((matriz_resultados[,k+1])-1)/N)
      yvetor_mc_resultados_sd[c]=sd(matriz_resultados[,k+1])
      yvetor_mc_resultados_media_pop[c]=mean(matriz_resultados[,k])
      yvetor_mc_resultados_media_est[c]=mean(matriz_resultados[,k+1])
    }
    
    
    
    ###################################################################
    ######################MEAN SQUARED ERROR###########################
    ###################################################################
    MSE_woold=WOOLxvetor_mc_resultados_vies^2+WOOLxvetor_mc_resultados_sd^2
    MSE_cic=xvetor_mc_resultados_vies^2+xvetor_mc_resultados_sd^2
    MSE_woold
    MSE_cic
    vetor_heatmap[cc]=mean(MSE_woold/MSE_cic)
    vetor_parametro1[cc]=escala
    vetor_parametro2[cc]=forma
    
    
    MSE_woold=WOOLyvetor_mc_resultados_vies^2+WOOLyvetor_mc_resultados_sd^2
    MSE_cic=yvetor_mc_resultados_vies^2+yvetor_mc_resultados_sd^2
    MSE_woold
    MSE_cic
    MSE_woold/MSE_cic
  }
}

#install.packages("reshape")
#install.packages("ggplot2")
library(reshape)
library(ggplot2)


dado_heatmap=data.frame(vetor_heatmap,as.factor(vetor_parametro1),as.factor(vetor_parametro2))
names(dado_heatmap)=c("MSE","escala","forma")

ggplot(dado_heatmap, aes(escala, forma, fill= MSE)) + 
  geom_tile() + ggtitle("(a) Com invari^ncia no tempo") +
  scale_fill_gradient(low="white", high="blue")
