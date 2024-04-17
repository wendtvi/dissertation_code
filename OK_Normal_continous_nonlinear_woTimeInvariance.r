#normal sem time inva
mc_function=function(N,parametro1,parametro2){
  n1=500
  TT=6
  q=4
  matriz_resultados=matrix(NA, ncol = (3*((TT-q+1)*2)),nrow = N)
  WOOLmatriz_resultados=matrix(NA, ncol = (3*((TT-q+1)*2)),nrow = N)
  CICmatriz_resultados=matrix(NA, ncol = (3*((TT-q+1)*2)),nrow = N)
  LB_CICmatriz_resultados=matrix(NA, ncol = (3*((TT-q+1)*2)),nrow = N)
  matriz_variancia_estimada4=matrix(NA, ncol = 3,nrow = N)
  matriz_variancia_estimada5=matrix(NA, ncol = 3,nrow = N)
  matriz_variancia_estimada6=matrix(NA, ncol = 3,nrow = N)
  
  for(p in 1:N){
    ##########################################################
    #####################CENÁRIO SIMULAÇÃO####################
    ##########################################################
    
    #Gerando covariável
    Z_cov_mean_tratados=vector()
    Z_cov_tratados=matrix(NA,ncol=TT,nrow=n1)
    for (k in 1:nrow(Z_cov_tratados)){
      Z_cov_tratados[k,]=rexp(TT,1)
      Z_cov_mean_tratados[k]=mean(Z_cov_tratados[k,])
    }
    
    Z_cov_mean_controle=vector()
    Z_cov_controle=matrix(NA,ncol=TT,nrow=n1)
    for (k in 1:nrow(Z_cov_controle)){
      Z_cov_controle[k,]=rexp(TT,1)
      Z_cov_mean_controle[k]=mean(Z_cov_controle[k,])
    }
    
    #Variável de tratamento
    V=rlogis(n1,0,1)
    D=vector()
    D=rep(-0.5,n1)+(Z_cov_mean_tratados-1)+V>0
    pD=sum(D)/length(D) #Incidência de tratamento na população
    D=as.numeric(D)
    
    ct_0_controle_pre=rep(0,n1)
    ct_0_controle_pos=rep(0,n1)
    ct_inf_controle_pre=rep(0,n1)
    ct_inf_controle_pos=rep(0,n1)
    ct_0_trat_pre=rep(0,n1)
    ct_0_trat_pos=rep(0,n1)
    ct_inf_trat_pre=rep(0,n1)
    ct_inf_trat_pos=rep(0,n1)
    
    #Modelo de diferenças em diferença para cada cohort (01 primero é grupo e segundo é tempo: neste caso é grupo de controle no período pós trat)
    f5=rep(1,n1)*D
    f6=rep(1,n1)*D
    matriz_estado_naotratamento=matrix(NA,nrow = n1,ncol = TT+1)
    matriz_estado_naotratamento[,TT+1]=t(D)
    for (k in 1:TT){
      for(i in 1:nrow(matriz_estado_naotratamento)){
        if (k<q && matriz_estado_naotratamento[i,TT+1]==0) matriz_estado_naotratamento[i,k]=as.numeric(-1+(Z_cov_controle[i,k]-1)*5/2-2*D[i]+(Z_cov_controle[i,k]-1)*D[i]/4+rnorm(1,rnorm(1,parametro1,1)+Z_cov_controle[i,k]*1/k,parametro2)+ct_0_controle_pre[i])/10
        if (k>=q&& matriz_estado_naotratamento[i,TT+1]==0) matriz_estado_naotratamento[i,k]=as.numeric(-1+(Z_cov_controle[i,k]-1)*5/2-2*D[i]+(Z_cov_controle[i,k]-1)*D[i]/4+rnorm(1,rnorm(1,parametro1,1)+Z_cov_controle[i,k]*1/k,parametro2)+ct_0_controle_pos[i])/10
        if (k<q && matriz_estado_naotratamento[i,TT+1]==1) matriz_estado_naotratamento[i,k]=as.numeric(-1+(Z_cov_tratados[i,k]-1)*5/2-2*D[i]+(Z_cov_tratados[i,k]-1)*D[i]/4+rnorm(1,rnorm(1,parametro1,1)+Z_cov_tratados[i,k]*1/k,parametro2)+ct_0_trat_pre[i])/10
        if (k>=q&& matriz_estado_naotratamento[i,TT+1]==1) matriz_estado_naotratamento[i,k]=as.numeric(-1+(Z_cov_tratados[i,k]-1)*5/2-2*D[i]+(Z_cov_tratados[i,k]-1)*D[i]/4+rnorm(1,rnorm(1,parametro1,1)+Z_cov_tratados[i,k]*1/k,parametro2)+ct_0_trat_pos[i])/10
        
      }
    }
    
    #Gerando matriz de dados latentes.
    matriz_X_estrela=matrix(NA, ncol = TT+1,nrow = n1)
    matriz_X_estrela[,TT+1]=t(D)
    for (k in 1:ncol(matriz_X_estrela)-1){
      for (i in 1:nrow(matriz_X_estrela)){
        if (k<q && matriz_X_estrela[i,TT+1]==0) matriz_X_estrela[i,k]=as.numeric(-1+(Z_cov_controle[i,k]-1)*5/2-2*D[i]+(Z_cov_controle[i,k]-1)*D[i]/4+rnorm(1,rnorm(1,parametro1,1)+Z_cov_controle[i,k]*1/k,parametro2)+ct_0_controle_pre[i])/10
        if (k>=q&& matriz_X_estrela[i,TT+1]==0) matriz_X_estrela[i,k]=as.numeric(-1+(Z_cov_controle[i,k]-1)*5/2-2*D[i]+(Z_cov_controle[i,k]-1)*D[i]/4+rnorm(1,rnorm(1,parametro1,1)+Z_cov_controle[i,k]*1/k,parametro2)+ct_0_controle_pos[i])/10
        if (k<q && matriz_X_estrela[i,TT+1]==1) matriz_X_estrela[i,k]=as.numeric(-1+(Z_cov_controle[i,k]-1)*5/2-2*D[i]+(Z_cov_controle[i,k]-1)*D[i]/4+rnorm(1,rnorm(1,parametro1,1)+Z_cov_tratados[i,k]*1/k,parametro2)+ct_0_trat_pre[i])/10
        if (k>=q&& matriz_X_estrela[i,TT+1]==1) matriz_X_estrela[i,k]=as.numeric(-1+(Z_cov_controle[i,k]-1)*5/2-2*D[i]+(Z_cov_controle[i,k]-1)*D[i]/4+rnorm(1,rnorm(1,parametro1,1)+Z_cov_tratados[i,k]*1/k,parametro2)+ct_0_trat_pos[i])/10
        if (k==4 && matriz_X_estrela[i,TT+1]==1)matriz_X_estrela[i,k]=as.numeric(-1+0.5+(Z_cov_tratados[i,k]-1)-2*D[i]+rnorm(1,rnorm(1,parametro1,1)+Z_cov_tratados[i,k]*1/k,parametro2)+ct_inf_trat_pos[i])/10
        if (k==5 && matriz_X_estrela[i,TT+1]==1)matriz_X_estrela[i,k]=as.numeric(-1+0.5+(Z_cov_tratados[i,k]-1)-2*D[i]+0.2*f5[i]+rnorm(1,rnorm(1,parametro1,1)+Z_cov_tratados[i,k]*1/k,parametro2)+ct_inf_trat_pos[i])/10
        if (k==6 && matriz_X_estrela[i,TT+1]==1)matriz_X_estrela[i,k]=as.numeric(-1+0.5+(Z_cov_tratados[i,k]-1)-2*D[i]+0.2*f5[i]+0.3*f6[i]+rnorm(1,rnorm(1,parametro1,1)+Z_cov_tratados[i,k]*1/k,parametro2)+ct_inf_trat_pos[i])/10
      }
    }
    
    ####################################################
    ##GERANDO VAR DO ESTADO NAO TRATADO#################
    ####################################################
    #GERANDO VARIAVEL X^*
    X_estrela_vetor_0_0=vector()
    X_estrela_vetor_1_0=vector()
    for (k in 1:TT){
      X_estrela_vetor_0_0[k]=mean((matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==0,k]))
      X_estrela_vetor_1_0[k]=mean((matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==1,k]))
    }
    
    #GERANDO VARIAVEL X
    
    X_vetor_0_0=matrix(NA,nrow = length(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==0,1]),ncol = TT)
    X_vetor_1_0=matrix(NA,nrow = length(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==1,1]),ncol = TT)
    for (k in 1:TT){
      X_vetor_0_0[,k]=exp(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==0,k])
      X_vetor_1_0[,k]=exp(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==1,k])
    }
    
    #GERANDO VARIAVEL Y
    Y_vetor_0_0=matrix(NA,nrow = n1*(1-pD),ncol = TT)
    Y_vetor_1_0=matrix(NA,nrow = n1*(pD),ncol = TT)
    for (k in 1:TT){
      Y_vetor_0_0[,k]=sum(exp(X_estrela_vetor_0_0[k]))
      Y_vetor_1_0[,k]=sum(exp(X_estrela_vetor_1_0[k]))
    }
    
    
    ####################################################
    ##GERANDO VAR DO ESTADO TRATADO#####################
    ####################################################
    #GERANDO VARIAVEL X^*
    X_estrela_vetor_0_inf=vector()
    X_estrela_vetor_1_inf=vector()
    for (k in 1:TT){
      X_estrela_vetor_0_inf[k]=(mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==0,k]))
      X_estrela_vetor_1_inf[k]=(mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,k]))
    }
    
    
    #GERANDO VARIAVEL X
    X_vetor_0_inf=matrix(NA,nrow = length(matriz_X_estrela[matriz_X_estrela[,TT+1]==0,1]),ncol = TT)
    X_vetor_1_inf=matrix(NA,nrow = length(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,1]),ncol = TT)
    for (k in 1:TT){
      X_vetor_0_inf[,k]=exp(matriz_X_estrela[matriz_X_estrela[,TT+1]==0,k])
      X_vetor_1_inf[,k]=exp(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,k])
    }
    
    #GERANDO VARIAVEL Y
    Y_vetor_0_inf=matrix(NA,nrow = n1*(1-pD),ncol = TT)
    Y_vetor_1_inf=matrix(NA,nrow = n1*(pD),ncol = TT)
    for (k in 1:TT){
      Y_vetor_0_inf[,k]=sum(exp(X_estrela_vetor_0_inf[k]))
      Y_vetor_1_inf[,k]=sum(exp(X_estrela_vetor_1_inf[k]))
    }
    
    ##########################################################
    #####################ESTIMADORES##########################
    ##########################################################
    
    G_Ybar_11_t4=log(mean(X_vetor_1_inf[,q]))
    G_Ybar_11_t5=log(mean(X_vetor_1_inf[,q+1]))
    G_Ybar_11_t6=log(mean(X_vetor_1_inf[,q+2]))
    
    G_Ybar_10_t3=log(mean(X_vetor_1_inf[,q-1]))
    G_Ybar_00_t3=log(mean(X_vetor_0_inf[,q-1]))
    G_Ybar_01_t4=log(mean(X_vetor_0_inf[,q]))
    G_Ybar_01_t5=log(mean(X_vetor_0_inf[,q+1]))
    G_Ybar_01_t6=log(mean(X_vetor_0_inf[,q+2]))
    
    
    ##########################################################
    #####################RESULTADOS###########################
    ##########################################################
    tau_4=X_estrela_vetor_1_inf[q]-X_estrela_vetor_1_0[q]
    tau_5=X_estrela_vetor_1_inf[q+1]-X_estrela_vetor_1_0[q+1]
    tau_6=X_estrela_vetor_1_inf[q+2]-X_estrela_vetor_1_0[q+2]
    
    tau_4_hat=X_estrela_vetor_1_inf[q]-mean(G_Ybar_01_t4-G_Ybar_00_t3+G_Ybar_10_t3)
    tau_5_hat=X_estrela_vetor_1_inf[q+1]-mean((G_Ybar_01_t5-G_Ybar_00_t3+G_Ybar_10_t3))
    tau_6_hat=X_estrela_vetor_1_inf[q+2]-mean((G_Ybar_01_t6-G_Ybar_00_t3+G_Ybar_10_t3))
    
    gamma_4=mean(X_vetor_1_inf[,q])-mean(X_vetor_1_0[,q])
    gamma_5=mean(X_vetor_1_inf[,q+1])-mean(X_vetor_1_0[,q+1])
    gamma_6=mean(X_vetor_1_inf[,q+2])-mean(X_vetor_1_0[,q+2])
    
    gamma4_hat=mean(X_vetor_1_inf[,q])-mean(exp((G_Ybar_01_t4-G_Ybar_00_t3+G_Ybar_10_t3)))
    gamma5_hat=mean(X_vetor_1_inf[,q+1])-mean(exp((G_Ybar_01_t5-G_Ybar_00_t3+G_Ybar_10_t3)))
    gamma6_hat=mean(X_vetor_1_inf[,q+2])-mean(exp((G_Ybar_01_t6-G_Ybar_00_t3+G_Ybar_10_t3)))
    
    
    kappa_4=Y_vetor_1_inf[q]-Y_vetor_1_0[q]
    kappa_5=Y_vetor_1_inf[q+1]-Y_vetor_1_0[q+1]
    kappa_6=Y_vetor_1_inf[q+2]-Y_vetor_1_0[q+2]
    
    
    kappa4_hat=Y_vetor_1_inf[q]-(log(G_Ybar_11_t4)*(n1*pD))
    kappa5_hat=Y_vetor_1_inf[q+1]-(log(G_Ybar_11_t5)*(n1*pD))
    kappa6_hat=Y_vetor_1_inf[q+2]-(log(G_Ybar_11_t6)*(n1*pD))
    
    WOOLmatriz_resultados[p,1]=tau_4
    WOOLmatriz_resultados[p,2]=tau_4_hat
    WOOLmatriz_resultados[p,3]=tau_5
    WOOLmatriz_resultados[p,4]=tau_5_hat
    WOOLmatriz_resultados[p,5]=tau_6
    WOOLmatriz_resultados[p,6]=tau_6_hat
    WOOLmatriz_resultados[p,7]=gamma_4
    WOOLmatriz_resultados[p,8]=gamma4_hat
    WOOLmatriz_resultados[p,9]=gamma_5
    WOOLmatriz_resultados[p,10]=gamma5_hat
    WOOLmatriz_resultados[p,11]=gamma_6
    WOOLmatriz_resultados[p,12]=gamma6_hat
    WOOLmatriz_resultados[p,13]=kappa_4
    WOOLmatriz_resultados[p,14]=kappa4_hat
    WOOLmatriz_resultados[p,15]=kappa_5
    WOOLmatriz_resultados[p,16]=kappa5_hat
    WOOLmatriz_resultados[p,17]=kappa_6
    WOOLmatriz_resultados[p,18]=kappa6_hat
    
    
    ##########################################################
    #####################ESTIMADORES##########################
    ##########################################################
    
    G_Ybar_11_t4=mean(log(X_vetor_1_inf[,q]))
    G_Ybar_11_t5=mean(log(X_vetor_1_inf[,q+1]))
    G_Ybar_11_t6=mean(log(X_vetor_1_inf[,q+2]))
    
    G_Ybar_10_t3=mean(log(X_vetor_1_inf[,q-1]))
    G_Ybar_00_t3=mean(log(X_vetor_0_inf[,q-1]))
    G_Ybar_01_t4=mean(log(X_vetor_0_inf[,q]))
    G_Ybar_01_t5=mean(log(X_vetor_0_inf[,q+1]))
    G_Ybar_01_t6=mean(log(X_vetor_0_inf[,q+2]))
    
    
    ##########################################################
    #####################RESULTADOS###########################
    ##########################################################
    tau_4=X_estrela_vetor_1_inf[q]-X_estrela_vetor_1_0[q]
    tau_5=X_estrela_vetor_1_inf[q+1]-X_estrela_vetor_1_0[q+1]
    tau_6=X_estrela_vetor_1_inf[q+2]-X_estrela_vetor_1_0[q+2]
    
    tau_4_hat=X_estrela_vetor_1_inf[q]-mean(G_Ybar_01_t4-G_Ybar_00_t3+G_Ybar_10_t3)
    tau_5_hat=X_estrela_vetor_1_inf[q+1]-mean((G_Ybar_01_t5-G_Ybar_00_t3+G_Ybar_10_t3))
    tau_6_hat=X_estrela_vetor_1_inf[q+2]-mean((G_Ybar_01_t6-G_Ybar_00_t3+G_Ybar_10_t3))
    
    gamma_4=exp(mean(log(X_vetor_1_inf[,q])))-exp(mean(log(X_vetor_1_0[,q])))
    gamma_5=exp(mean(log(X_vetor_1_inf[,q+1])))-exp(mean(log(X_vetor_1_0[,q+1])))
    gamma_6=exp(mean(log(X_vetor_1_inf[,q+2])))-exp(mean(log(X_vetor_1_0[,q+2])))
    
    gamma4_hat=exp(mean(log(X_vetor_1_inf[,q])))-(exp((G_Ybar_01_t4-G_Ybar_00_t3+G_Ybar_10_t3)))
    gamma5_hat=exp(mean(log(X_vetor_1_inf[,q+1])))-(exp((G_Ybar_01_t5-G_Ybar_00_t3+G_Ybar_10_t3)))
    gamma6_hat=exp(mean(log(X_vetor_1_inf[,q+2])))-(exp((G_Ybar_01_t6-G_Ybar_00_t3+G_Ybar_10_t3)))
    
    
    kappa_4=Y_vetor_1_inf[q]-Y_vetor_1_0[q]
    kappa_5=Y_vetor_1_inf[q+1]-Y_vetor_1_0[q+1]
    kappa_6=Y_vetor_1_inf[q+2]-Y_vetor_1_0[q+2]
    
    
    kappa4_hat=Y_vetor_1_inf[q]-(log(G_Ybar_11_t4)*(n1*pD))
    kappa5_hat=Y_vetor_1_inf[q+1]-(log(G_Ybar_11_t5)*(n1*pD))
    kappa6_hat=Y_vetor_1_inf[q+2]-(log(G_Ybar_11_t6)*(n1*pD))
    
    WOOLmatriz_resultados[p,1]=tau_4
    WOOLmatriz_resultados[p,2]=tau_4_hat
    WOOLmatriz_resultados[p,3]=tau_5
    WOOLmatriz_resultados[p,4]=tau_5_hat
    WOOLmatriz_resultados[p,5]=tau_6
    WOOLmatriz_resultados[p,6]=tau_6_hat
    WOOLmatriz_resultados[p,7]=gamma_4
    WOOLmatriz_resultados[p,8]=gamma4_hat
    WOOLmatriz_resultados[p,9]=gamma_5
    WOOLmatriz_resultados[p,10]=gamma5_hat
    WOOLmatriz_resultados[p,11]=gamma_6
    WOOLmatriz_resultados[p,12]=gamma6_hat
    WOOLmatriz_resultados[p,13]=kappa_4
    WOOLmatriz_resultados[p,14]=kappa4_hat
    WOOLmatriz_resultados[p,15]=kappa_5
    WOOLmatriz_resultados[p,16]=kappa5_hat
    WOOLmatriz_resultados[p,17]=kappa_6
    WOOLmatriz_resultados[p,18]=kappa6_hat
    
    
    ##########################################################
    #####################ESTIMADORES##########################
    ##########################################################
    Y10=exp((matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q-1]))#variável resposta observada para grupo dos tratados no período pré tratamento t=3
    F_Y10=plnorm((Y10), meanlog = (mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==0,q-1])),sdlog=sd(matriz_X_estrela[matriz_X_estrela[,TT+1]==0,q-1]))
    
    F_inver_F_Y10_t4=qlnorm(F_Y10, meanlog = (mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==0,q])),sdlog=sd(matriz_X_estrela[matriz_X_estrela[,TT+1]==0,q]))
    F_inver_F_Y10_t5=qlnorm(F_Y10, meanlog = (mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==0,q+1])),sdlog=sd(matriz_X_estrela[matriz_X_estrela[,TT+1]==0,q+1]))
    F_inver_F_Y10_t6=qlnorm(F_Y10, meanlog = (mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==0,q+2])),sdlog=sd(matriz_X_estrela[matriz_X_estrela[,TT+1]==0,q+2]))
    
    F_inver_F_Y10_t4=mean(log(F_inver_F_Y10_t4))
    F_inver_F_Y10_t5=mean(log(F_inver_F_Y10_t5))
    F_inver_F_Y10_t6=mean(log(F_inver_F_Y10_t6))
    
    
    ##########################################################
    #####################RESULTADOS###########################
    ##########################################################
    tau_4=X_estrela_vetor_1_inf[q]-X_estrela_vetor_1_0[q]
    tau_5=X_estrela_vetor_1_inf[q+1]-X_estrela_vetor_1_0[q+1]
    tau_6=X_estrela_vetor_1_inf[q+2]-X_estrela_vetor_1_0[q+2]
    
    tau_4_hat=X_estrela_vetor_1_inf[q]-((F_inver_F_Y10_t4))
    tau_5_hat=X_estrela_vetor_1_inf[q+1]-((F_inver_F_Y10_t5))
    tau_6_hat=X_estrela_vetor_1_inf[q+2]-((F_inver_F_Y10_t6))
    
    gamma_4=exp(mean(log(X_vetor_1_inf[,q])))-exp(mean(log(X_vetor_1_0[,q])))
    gamma_5=exp(mean(log(X_vetor_1_inf[,q+1])))-exp(mean(log(X_vetor_1_0[,q+1])))
    gamma_6=exp(mean(log(X_vetor_1_inf[,q+2])))-exp(mean(log(X_vetor_1_0[,q+2])))
    
    gamma4_hat=exp(mean(log(X_vetor_1_inf[,q])))-exp(F_inver_F_Y10_t4)
    gamma5_hat=exp(mean(log(X_vetor_1_inf[,q+1])))-exp(F_inver_F_Y10_t5)
    gamma6_hat=exp(mean(log(X_vetor_1_inf[,q+2])))-exp(F_inver_F_Y10_t6)
    
    kappa_4=Y_vetor_1_inf[q]-Y_vetor_1_0[q]
    kappa_5=Y_vetor_1_inf[q+1]-Y_vetor_1_0[q+1]
    kappa_6=Y_vetor_1_inf[q+2]-Y_vetor_1_0[q+2]
    
    
    kappa4_hat=Y_vetor_1_inf[q]-F_inver_F_Y10_t4
    kappa5_hat=Y_vetor_1_inf[q+1]-F_inver_F_Y10_t5
    kappa6_hat=Y_vetor_1_inf[q+2]-F_inver_F_Y10_t6
    
    
    
    matriz_resultados[p,1]=tau_4
    matriz_resultados[p,2]=tau_4_hat
    matriz_resultados[p,3]=tau_5
    matriz_resultados[p,4]=tau_5_hat
    matriz_resultados[p,5]=tau_6
    matriz_resultados[p,6]=tau_6_hat
    matriz_resultados[p,7]=gamma_4
    matriz_resultados[p,8]=gamma4_hat
    matriz_resultados[p,9]=gamma_5
    matriz_resultados[p,10]=gamma5_hat
    matriz_resultados[p,11]=gamma_6
    matriz_resultados[p,12]=gamma6_hat
    matriz_resultados[p,13]=kappa_4
    matriz_resultados[p,14]=kappa4_hat
    matriz_resultados[p,15]=kappa_5
    matriz_resultados[p,16]=kappa5_hat
    matriz_resultados[p,17]=kappa_6
    matriz_resultados[p,18]=kappa6_hat
    
    
    
    
    ##########################################################
    #####################ESTIMADORES CIC######################
    ##########################################################
    Y_10=exp((matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q-1]))#variável resposta observada para grupo dos tratados no período pré tratamento t=3
    Y_11_t4=exp((matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q]))
    Y_11_t5=exp((matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q+1]))
    Y_11_t6=exp((matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q+2]))
    Y_00=exp((matriz_X_estrela[matriz_X_estrela[,TT+1]==0,q-1]))
    Y_01_t4=exp((matriz_X_estrela[matriz_X_estrela[,TT+1]==0,q]))
    Y_01_t5=exp((matriz_X_estrela[matriz_X_estrela[,TT+1]==0,q+1]))
    Y_01_t6=exp((matriz_X_estrela[matriz_X_estrela[,TT+1]==0,q+2]))
    
    #UB
    F00_hat_LB=vector()
    for (i in 1:length(Y_10)) { F00_hat_LB[i]=mean(Y_00<=Y_10[i]) }
    F_inver_01_hat_LB_t4=vector()
    F_inver_01_hat_LB_t5=vector()
    F_inver_01_hat_LB_t6=vector()
    F01_hat_LB_t4=vector()
    for (i in 1:length(Y_01_t4)) { F01_hat_LB_t4[i]=mean(Y_01_t4<=Y_01_t4[i]) }
    for (i in 1:length(F00_hat_LB)) { F_inver_01_hat_LB_t4[i]=min(Y_01_t4[F01_hat_LB_t4>=F00_hat_LB[i]]) }
    F01_hat_LB_t5=vector()
    for (i in 1:length(Y_01_t5)) { F01_hat_LB_t5[i]=mean(Y_01_t5<=Y_01_t5[i]) }
    for (i in 1:length(F00_hat_LB)) { F_inver_01_hat_LB_t5[i]=min(Y_01_t5[F01_hat_LB_t5>=F00_hat_LB[i]]) }
    F01_hat_LB_t6=vector()
    for (i in 1:length(Y_01_t6)) { F01_hat_LB_t6[i]=mean(Y_01_t6<=Y_01_t6[i]) }
    for (i in 1:length(F00_hat_LB)) { F_inver_01_hat_LB_t6[i]=min(Y_01_t6[F01_hat_LB_t6>=F00_hat_LB[i]]) }
    ATT_CIC_LB_t4=exp(mean(log(Y_11_t4)))-exp(mean(log(F_inver_01_hat_LB_t4)))
    ATT_CIC_LB_t5=exp(mean(log(Y_11_t5)))-exp(mean(log(F_inver_01_hat_LB_t5)))
    ATT_CIC_LB_t6=exp(mean(log(Y_11_t6)))-exp(mean(log(F_inver_01_hat_LB_t6)))
    
    
    
    ##########################################################
    #####################RESULTADOS CIC#######################
    ##########################################################
    tau_4=X_estrela_vetor_1_inf[q]-X_estrela_vetor_1_0[q]
    tau_5=X_estrela_vetor_1_inf[q+1]-X_estrela_vetor_1_0[q+1]
    tau_6=X_estrela_vetor_1_inf[q+2]-X_estrela_vetor_1_0[q+2]
    
    tau_4_hat=X_estrela_vetor_1_inf[q]-((F_inver_F_Y10_t4))
    tau_5_hat=X_estrela_vetor_1_inf[q+1]-((F_inver_F_Y10_t5))
    tau_6_hat=X_estrela_vetor_1_inf[q+2]-((F_inver_F_Y10_t6))
    
    gamma_4=exp(mean(log(X_vetor_1_inf[,q])))-exp(mean(log(X_vetor_1_0[,q])))
    gamma_5=exp(mean(log(X_vetor_1_inf[,q+1])))-exp(mean(log(X_vetor_1_0[,q+1])))
    gamma_6=exp(mean(log(X_vetor_1_inf[,q+2])))-exp(mean(log(X_vetor_1_0[,q+2])))
    
    gamma4_hat=ATT_CIC_LB_t4
    gamma5_hat=ATT_CIC_LB_t5
    gamma6_hat=ATT_CIC_LB_t6
    
    kappa_4=Y_vetor_1_inf[q]-Y_vetor_1_0[q]
    kappa_5=Y_vetor_1_inf[q+1]-Y_vetor_1_0[q+1]
    kappa_6=Y_vetor_1_inf[q+2]-Y_vetor_1_0[q+2]
    
    
    kappa4_hat=Y_vetor_1_inf[q]-F_inver_F_Y10_t4
    kappa5_hat=Y_vetor_1_inf[q+1]-F_inver_F_Y10_t5
    kappa6_hat=Y_vetor_1_inf[q+2]-F_inver_F_Y10_t6
    
    
    
    CICmatriz_resultados[p,1]=tau_4
    CICmatriz_resultados[p,2]=tau_4_hat
    CICmatriz_resultados[p,3]=tau_5
    CICmatriz_resultados[p,4]=tau_5_hat
    CICmatriz_resultados[p,5]=tau_6
    CICmatriz_resultados[p,6]=tau_6_hat
    CICmatriz_resultados[p,7]=gamma_4
    CICmatriz_resultados[p,8]=gamma4_hat
    CICmatriz_resultados[p,9]=gamma_5
    CICmatriz_resultados[p,10]=gamma5_hat
    CICmatriz_resultados[p,11]=gamma_6
    CICmatriz_resultados[p,12]=gamma6_hat
    CICmatriz_resultados[p,13]=kappa_4
    CICmatriz_resultados[p,14]=kappa4_hat
    CICmatriz_resultados[p,15]=kappa_5
    CICmatriz_resultados[p,16]=kappa5_hat
    CICmatriz_resultados[p,17]=kappa_6
    CICmatriz_resultados[p,18]=kappa6_hat
    
    
  }
  return(list(WOOLmatriz_resultados,matriz_resultados,CICmatriz_resultados,LB_CICmatriz_resultados,matriz_variancia_estimada4,matriz_variancia_estimada5,matriz_variancia_estimada6))
}
