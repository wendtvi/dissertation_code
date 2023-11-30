par(mfrow=c(1,2))

######################################
###########IMAGEM 1 DISSERTAÇAO#######
######################################

#Graf tendencia paralela nos dados (continuo)
t=seq(1,10,.1)
y_00=rlnorm(n = length(t), meanlog = 0.5, sdlog = .01)
y_01_t=rlnorm(n = length(t), meanlog = (0.5+2*t), sdlog = .01)
y_10=rlnorm(n = length(t), meanlog = (0.5+2), sdlog = .01)
y_11_t=rlnorm(n = length(t), meanlog = (0.5+2*t+2), sdlog = .001)
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


