require(ggplot2)

# PARAMETERS
g<-function(q) (1-q)^4

optmfun<-function(N) {

  # ACTION FUNCTION
  A<-function(x1,x2) {
    if (x1==x2) {
      x2
    } else {
      integrate(function(q)q*g(q),x1,x2)$value/integrate(g,x1,x2)$value
    }
  }
  X<-function(x1,x2) {
    if (x1==x2) {
      x1
    } else {
      uniroot(function(x3)x2-(A(x1,x2)+A(x2,x3))/2,c(x2,N+1))$root
    }
  }
  
  # NONLINEAR SHOOTING
  shooter<-function(x0,x1) {
    x<-rep(0,N+2)
    x[1]<-x0
    x[2]<-x1
    for (j in 1:N) {
      x[j+2]<-X(x[j],x[j+1])
    }
    return(x)
  }	
  
  # COMPUTE KNOTS
  x1.sol<-uniroot(function(x1)shooter(0,x1)[N+2]-1,c(0,1))
  print(N)
  print(x1.sol)
  sol<-data.frame(x.sol=shooter(0,x1.sol$root),y.sol=c(seq(0,1,1/N),1))
}

sol0<-optmfun(100)

# PLOT 
ggplot()+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
  geom_step(aes(x=sol0$y.sol,y=sol0$x.sol,colour="red"))+
  geom_line(aes(x=sol0$x.sol^2,y=sol0$x.sol,colour="black"))+
  labs(xlab("m_{k}"),ylab("q_{k}"))
  theme_bw()
ggsave("kozelplot.pdf",width=5,height=5)
