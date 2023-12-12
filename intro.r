library(latex2exp)
options(scipen = 999)

par(mfrow=c(1,2))

######################################
###########IMAGEM 1 DISSERTAÇAO#######
######################################

#Graf tendencia paralela nos dados (continuo)
t=seq(1,10,.1)
y_00=rlnorm(n = length(t), meanlog = 0.5, sdlog = 1)
y_01_t=rlnorm(n = length(t), meanlog = (0.5+2*t), sdlog = 1)
y_10=rlnorm(n = length(t), meanlog = (0.5+2), sdlog = 1+t)
y_11_t=rlnorm(n = length(t), meanlog = (0.5+2*t+2), sdlog = 1+t)
T_star=10

(y_01_t-y_00+y_10)
plot(x=t,y=(y_01_t-y_00+y_10),type = 'l',ylab = "")
lines(y=((y_11_t)),x=t,col="red")
lines(y=((y_11_t)),x=rep(6, length(y_11_t)),col="gray")
legend('topleft', legend=c(TeX('$Y_{i,11}$'),TeX('contrafactual')), 
       fill=c("red","black"), ncol=1, xpd=NA)

(log(y_01_t)-log(y_00)+log(y_10))
plot(x=t,y=(log(y_01_t)-log(y_00)+log(y_10)),type = 'l',ylab = "")
lines(y=(log(y_11_t)),x=t,col="red")
lines(y=(log(y_11_t)),x=rep(6, length(y_11_t)),col="gray")
legend('topleft', legend=c(TeX('$Z_{i,11}$'),TeX('contrafactual')), 
       fill=c("red","black"), ncol=1, xpd=NA)



Qy_00=ecdf((y_00-mean(y_00))/sd(y_00))
Qy_01_t=ecdf((y_01_t-mean(y_01_t))/sd(y_01_t))
Qy_10=ecdf((y_10-mean(y_10))/sd(y_10))
Qy_11_t=ecdf((y_11_t-mean(y_11_t))/sd(y_11_t))

par(mfrow=c(2,2))
plot(Qy_00,ylab = TeX('$F_{Y_{00}}$'),main="",xlab="Y")
plot(Qy_01_t,ylab = TeX('$F_{Y_{01}}$'),main="",xlab="Y")
plot(Qy_10,ylab = TeX('$F_{Y_{10}}$'),main="",xlab="Y")
plot(Qy_11_t,ylab = TeX('$F_{Y_{11}}$'),main="",xlab="Y")

#Supondo distribuição conhecida

quantile_y00=sort((y_00-mean(y_00))/sd(y_00))
quantile_y01=sort((y_01_t-mean(y_01_t))/sd(y_01_t))
quantile_y10=sort((y_10-mean(y_10))/sd(y_10))
quantile_y11=sort((y_11_t-mean(y_11_t))/sd(y_11_t))

e_cdf00 = plnorm(y_10,meanlog = mean(log(y_00)),sdlog = sd(log(y_00)))
e_cdf01 =  plnorm(y_10,meanlog = mean(log(y_01_t)),sdlog = sd(log(y_01_t)))
e_cdf10 =  plnorm(y_10,meanlog = mean(log(y_10)),sdlog = sd(log(y_10)))
e_cdf11 =  plnorm(y_10,meanlog = mean(log(y_11_t)),sdlog = sd(log(y_11_t)))

q_00=qlnorm(e_cdf10,meanlog = mean(log(y_00)),sdlog = sd(log(y_00)))
q_01=qlnorm(e_cdf10,meanlog = mean(log(y_01_t)),sdlog = sd(log(y_01_t)))
q_10=qlnorm(e_cdf10,meanlog = mean(log(y_10)),sdlog = sd(log(y_10)))
q_11=qlnorm(e_cdf10,meanlog = mean(log(y_11_t)),sdlog = sd(log(y_11_t)))

par(mfrow=c(2,1))
#tendencia paralela quantilica
#Dist conhecida
plot(abs((q_00-mean(q_00))/sd(q_00)-(q_01-mean(q_01))/sd(q_01))-abs((q_10-mean(q_10))/sd(q_10)-(q_11-mean(q_11))/sd(q_11)),,type = 'l',ylim = c(-2,2))
mean(abs((q_00-mean(q_00))/sd(q_00)-(q_01-mean(q_01))/sd(q_01))-abs((q_10-mean(q_10))/sd(q_10)-(q_11-mean(q_11))/sd(q_11)))

#roth e santanna
e_cdf00 = plnorm(y_10,meanlog = mean(log(y_00)),sdlog = sd(log(y_00)))
e_cdf01 =  plnorm(y_10,meanlog = mean(log(y_01_t)),sdlog = sd(log(y_01_t)))
e_cdf10 =  plnorm(y_10,meanlog = mean(log(y_10)),sdlog = sd(log(y_10)))
e_cdf11 =  plnorm(y_10,meanlog = mean(log(y_11_t)),sdlog = sd(log(y_11_t)))
plot(abs((e_cdf00-mean(e_cdf00))/sd(e_cdf00)-(e_cdf01-mean(e_cdf01))/sd(e_cdf01))-abs((e_cdf10-mean(e_cdf10))/sd(e_cdf10)-(e_cdf11-mean(e_cdf11))/sd(e_cdf11)),,type = 'l',ylim = c(-2,2))
mean(abs((e_cdf00-mean(e_cdf00))/sd(e_cdf00)-(e_cdf01-mean(e_cdf01))/sd(e_cdf01))-abs((e_cdf10-mean(e_cdf10))/sd(e_cdf10)-(e_cdf11-mean(e_cdf11))/sd(e_cdf11)))

#tendencia paralela quantilica
#Distribuição empírica
fun_ecdf=ecdf(x)
Py_00=fun_ecdf((y_00-mean(y_00))/sd(y_00))
Py_01_t=fun_ecdf((y_01_t-mean(y_01_t))/sd(y_01_t))
Py_10=fun_ecdf((y_10-mean(y_10))/sd(y_10))
Py_11_t=fun_ecdf((y_11_t-mean(y_11_t))/sd(y_11_t))

quantile_y00=sort((y_00-mean(y_00))/sd(y_00))
quantile_y01=sort((y_01_t-mean(y_01_t))/sd(y_01_t))
quantile_y10=sort((y_10-mean(y_10))/sd(y_10))
quantile_y11=sort((y_11_t-mean(y_11_t))/sd(y_11_t))

e_cdf00 <- 1:length(quantile_y00) / length(quantile_y00)
e_cdf01 <- 1:length(quantile_y01) / length(quantile_y01)
e_cdf10 <- 1:length(quantile_y10) / length(quantile_y10)
e_cdf11 <- 1:length(quantile_y11) / length(quantile_y11)

prob=seq(0,1,0.01)
c=0
for(p in prob){
  c=c+1
  diff_quantiles0[c]=abs(quantile_y00[which(e_cdf00  >= p)[1]]-quantile_y01[which(e_cdf01 >= p)[1]])
  diff_quantiles1[c]=abs(quantile_y10[which(e_cdf10 >= p)[1]]-quantile_y11[which(e_cdf11 >= p)[1]])
  
}

plot(diff_quantiles0-diff_quantiles1,type = 'l',ylim = c(-2,2),x=prob,ylab="Diferença entre quantis")
mean(diff_quantiles0-diff_quantiles1)














#Graf tendencia paralela nos dados (discreto)
t=seq(1,10,.1)
y_00=rpois(n = length(t),lambda =  exp(0.5))
y_01_t=rpois(n = length(t), lambda = exp(0.5+2*t))
y_10=rpois(n = length(t), lambda = exp(0.5+2))
y_11_t=rpois(n = length(t), lambda = exp(0.5+2*t+2))
T_star=10

(y_01_t-y_00+y_10)
plot(x=t,y=(y_01_t-y_00+y_10),type = 'l',ylab = "")
lines(y=((y_11_t)),x=t,col="red")
lines(y=((y_11_t)),x=rep(6, length(y_11_t)),col="gray")
legend('topleft', legend=c(TeX('$Y_{i,11}$'),TeX('contrafactual')), 
       fill=c("red","black"), ncol=1, xpd=NA)

(log(y_01_t+1)-log(y_00+1)+log(y_10+1))
plot(x=t,y=(log(y_01_t+1)-log(y_00+1)+log(y_10+1)),type = 'l',ylab = "")
lines(y=(log(y_11_t+1)),x=t,col="red")
lines(y=(log(y_11_t+1)),x=rep(6, length(y_11_t)),col="gray")
legend('topleft', legend=c(TeX('$Z_{i,11}$'),TeX('contrafactual')), 
       fill=c("red","black"), ncol=1, xpd=NA)


Qy_00=ecdf((y_00-mean(y_00))/sd(y_00))
Qy_01_t=ecdf((y_01_t-mean(y_01_t))/sd(y_01_t))
Qy_10=ecdf((y_10-mean(y_10))/sd(y_10))
Qy_11_t=ecdf((y_11_t-mean(y_11_t))/sd(y_11_t))

par(mfrow=c(2,2))
plot(Qy_00,ylab = TeX('$F_{Y_{00}}$'),main="",xlab="Y")
plot(Qy_01_t,ylab = TeX('$F_{Y_{01}}$'),main="",xlab="Y")
plot(Qy_10,ylab = TeX('$F_{Y_{10}}$'),main="",xlab="Y")
plot(Qy_11_t,ylab = TeX('$F_{Y_{11}}$'),main="",xlab="Y")



fun_ecdf=ecdf(x)
Qy_00=fun_ecdf((y_00-mean(y_00))/sd(y_00))
Qy_01_t=fun_ecdf((y_01_t-mean(y_01_t))/sd(y_01_t))
Qy_10=fun_ecdf((y_10-mean(y_10))/sd(y_10))
Qy_11_t=fun_ecdf((y_11_t-mean(y_11_t))/sd(y_11_t))


quantile_y00=sort((y_00-mean(y_00))/sd(y_00))
quantile_y01=sort((y_01_t-mean(y_01_t))/sd(y_01_t))
quantile_y10=sort((y_10-mean(y_10))/sd(y_10))
quantile_y11=sort((y_11_t-mean(y_11_t))/sd(y_11_t))

e_cdf00 <- 1:length(quantile_y00) / length(quantile_y00)
e_cdf01 <- 1:length(quantile_y01) / length(quantile_y01)
e_cdf10 <- 1:length(quantile_y10) / length(quantile_y10)
e_cdf11 <- 1:length(quantile_y11) / length(quantile_y11)

prob=seq(0,1,0.01)
c=0
for(p in prob){
  c=c+1
  diff_quantiles0[c]=quantile_y00[which(e_cdf00  >= p)[1]]-quantile_y01[which(e_cdf01 >= p)[1]]
  diff_quantiles1[c]=quantile_y10[which(e_cdf10 >= p)[1]]-quantile_y11[which(e_cdf11 >= p)[1]]
  
}

par(mfrow=c(1,1))

plot(diff_quantiles0-diff_quantiles1,type = 'l',ylim = c(-2,2),x=prob)

