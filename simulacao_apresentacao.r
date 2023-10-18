par(mfrow=c(1,2))

#Roth e santanna problema
x0=qnorm(seq(0,1,0.01),0,1)
kx0=qnorm(seq(0,1,0.01),1,1.0)

x1=qnorm(seq(0,1,0.01),0.5,1)
kx1=qnorm(seq(0,1,0.01),1.5,1.5)

plot(y=seq(0,1,0.01),x=x0,type = 'l',xlim = c(-3,3))
lines(y=seq(0,1,0.01),x=kx0,col=10)
lines(y=seq(0,1,0.01),x=x1,col=11)
lines(y=seq(0,1,0.01),x=kx1,col=15)


plot(y=seq(0,1,0.01),x=x1,type = 'l',xlim = c(-3,3))
lines(y=seq(0,1,0.01),x=kx1)


#Minha ideia
x0=qnorm(seq(0,1,0.01),0,1)
kx0=qnorm(seq(0,1,0.01),1,1.0)

x1=qnorm(seq(0,1,0.01),0.5,1)
kx1=qnorm(seq(0,1,0.01),1.5,1.5)

plot(y=seq(0,1,0.01),x=x0,type = 'l',xlim = c(-3,3))
lines(y=seq(0,1,0.01),x=kx0,col=10)
lines(y=seq(0,1,0.01),x=x1,col=11)
lines(y=seq(0,1,0.01),x=kx1,col=15)


plot(y=seq(0,1,0.01),x=x1,type = 'l',xlim = c(-3,3))
lines(y=seq(0,1,0.01),x=kx1)
