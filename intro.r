library(latex2exp)

par(mfrow=c(1,2))

######################################
###########IMAGEM 1 DISSERTAÇAO#######
######################################

#Graf tendencia paralela nos dados (continuo)
t=seq(1,10,.1)
y_00=rlnorm(n = length(t), meanlog = 0.5, sdlog = .01)
y_01_t=rlnorm(n = length(t), meanlog = (0.5+2*t), sdlog = .01)
y_10=rlnorm(n = length(t), meanlog = (0.5+2), sdlog = .01)
y_11_t=rlnorm(n = length(t), meanlog = (0.5+2*t+2), sdlog = .01)
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



fun_ecdf=ecdf(x)
Py_00=fun_ecdf((y_00-mean(y_00))/sd(y_00))
Py_01_t=fun_ecdf((y_01_t-mean(y_01_t))/sd(y_01_t))
Py_10=fun_ecdf((y_10-mean(y_10))/sd(y_10))
Py_11_t=fun_ecdf((y_11_t-mean(y_11_t))/sd(y_11_t))


Qy_00=(y_00-mean(y_00))/sd(y_00)[fun_ecdf((y_00-mean(y_00))/sd(y_00))]
Qy_01_t=fun_ecdf((y_01_t-mean(y_01_t))/sd(y_01_t))
Qy_10=fun_ecdf((y_10-mean(y_10))/sd(y_10))
Qy_11_t=fun_ecdf((y_11_t-mean(y_11_t))/sd(y_11_t))

par(mfrow=c(1,1))
plot((Qy_01_t-Qy_00)-(Qy_11_t-Qy_10),x=t,type = "l",ylim = c(-1,1))




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

par(mfrow=c(1,1))
plot((Qy_01_t-Qy_00)-(Qy_11_t-Qy_10),x=t,type = "l",
     ylab =TeX('($F_{Y_{00}}-F_{Y_{01}})-($F_{Y_{11}}-F_{Y_{10}})$'))

