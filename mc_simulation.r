##############################################
###############MCMC###########################
##############################################
#install.packages("devtools")
library(devtools)

#source_url("https://github.com/wendtvi/dissertation_code/blob/main/continous_nonlinear.r?raw=TRUE")

matriz_resultados=mc_function(N=1000)
vetor_mc_resultados_vies=vector()
vetor_mc_resultados_sd=vector()
vetor_mc_resultados_media_pop=vector()
vetor_mc_resultados_media_est=vector()
c=0
for (k in seq(1,ncol(matriz_resultados),2)){
  c=c+1
  vetor_mc_resultados_vies[c]=mean(matriz_resultados[,k])-mean(matriz_resultados[,k+1])
  vetor_mc_resultados_sd[c]=sd(matriz_resultados[,k+1])
  vetor_mc_resultados_media_pop[c]=mean(matriz_resultados[,k])
  vetor_mc_resultados_media_est[c]=mean(matriz_resultados[,k+1])
}

par(mfrow=c(1,2))
hist(matriz_resultados[,4],main="Distribuição estimador Wooldridge t=4 e sd 2 no modelo DiD",xlim=c(-2,2))
hist(matriz_resultados[,4],main="Distribuição estimador CIC t=4 e sd 2 no modelo DiD",xlim=c(-2,2))
