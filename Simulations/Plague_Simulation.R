 
rm(list=ls())

library(deSolve)
library(nls2)    # pas obligatoire (on peut utiliser nls au lieu de nls2)
library(sensitivity)
# ------------------ DONNEES

# modélisation par compartiments : 

# V=vecteur d'états, p=vecteur de paramètres


signal <-function(t){          # Fonction de température 
#  res<- 4*(cos(2*pi*t/15)+(1/9)*cos(6*pi*t/15)+(1/25)*cos(10*pi*t/15))  # les trois premiers termes d'un signal triangulaire 
# res<- (20/pi)*(sin(2*pi*t/12)+0.33*sin(6*pi*t/12)+0.2*sin(10*pi*t/12))
  res<-10*cos(2*pi*t/12)
  return(res)        
  
}


deriv <- function(t,y,p){
  with (as.list(c(p, y)), {
    P=S+I+R
    dS <- ((r*P)/(1+P/K)-m*S-c*(1-omega)*pi0*(0.75-0.25*tanh(teta+signal(t)-80))*(S/P)*F+epsilon*m2*I)   
    dI <- (c*(1-omega)*pi0*(0.75-0.25*tanh(teta+signal(t)-80))*(S/P)*F-m2*I )
    dR <- (epsilon2*m2*I-m*R)
    dF <- (f*pi2*(1-epsilon-epsilon2)*m2*I-c*F)
    dH <- (c*omega*pi0*(0.75-0.25*tanh(teta+signal(t)-80))*F-a*H)
    #cat("\n t=", t," dS=", dS, " dI=", dI," dR"=dR," dH"=dH)
    list(c(dS,dI,dR,dF,dH)) 
    
    
  })
}
 #  SIMULATION
# test d'une simulation :
N = 500                     # nombre d'élèves
Tmax = 120               # temps maximal
par <- list(r=0.4,K=50000,m=0.03, c=30,omega=0.02, pi0=0.9, teta=80,m2=30,epsilon=0.1,epsilon2=0.1,f=4,pi2=0.9,a=4,sigma=0.9) # vecteurs de paramètres 
V <- c(S=N-1, I=1,R=0,F=100 ,H=0)              # valeurs initiales du vecteur d'états
comparts <- ode(V, seq(1,Tmax,by=0.1), deriv, parms=par)
par(mfrow=c(2,2),mar=c(4,5,0.5,1))
plot(comparts[,"time"], comparts[,"H"],type="l",xlab="Temps",ylab="H(t)", lwd=2, col="red")
plot(comparts[,"time"], comparts[,"I"],type="l",xlab="Temps",ylab="I(t)", lwd=2, col="blue")
plot(comparts[,"time"], comparts[,"F"],type="l",xlab="Temps",ylab="F(t)", lwd=2, col="brown")
plot(comparts[,"time"], comparts[,"R"],type="l",xlab="Temps",ylab="R(t)", lwd=2, col="black")
#Bruit <- rnorm(Tmax,mean=0,sd=0.5)  #Bruit normal centré
Bruit <- sapply(rnorm(length(seq(1,Tmax,by=0.1)),mean=0,sd=1), function(x) max(x,0))  #Bruit positif (On prend le max entre une bruit normal et 0)


d <- data.frame(time=seq(1,Tmax,by=0.1),H=(comparts[,"H"]+Bruit))
save(d,file="NewData.Rda")
write.table(d, "tot.txt")


#  ESTIMATION
# Fonction pour l'estimation des paramètres
# Retourne un vecteur I simulé, qui doit être de même taille que celui des données.
model_tofit <- function(r,K,m,c,omega,pi0,teta,m2,epsilon,epsilon2,f,pi2,a,sigma, times){
  V <- c(S=N-1, I=1,R=0,F=100 ,H=0)
  #cat(r,K,m,c,omega,pi0,teta,m2,epsilon,epsilon2,f,pi2,a,sigma, "\n")
  Fsim <- ode(V, times, deriv, parms=list("r"=r,"K"=K,"m"=m,"c"=c,"omega"=omega, "pi0"=pi0,"teta"=teta,"m2"=m2,"epsilon"=epsilon,"epsilon2"=epsilon2,"f"=f,"pi2"=pi2,"a"=a,"sigma"=sigma))
  return(Fsim[,"H"])  
}

# estimation des paramètres p :                        
obs <- d$H 
res<-nls2(obs~model_tofit(r,K,m,c,omega,pi0,teta,m2,epsilon,epsilon2,f,pi2,a,sigma, T), start=list(pi0=0.2,a=2), data=list(T=d$time,r=0.4,K=50000, c=30,omega=0.9, teta=80,m2=30,epsilon=0.1,epsilon2=0.1,m=0.03,pi2=0.9,f=4,sigma=0.9), trace=TRUE, nls.control(warnOnly=TRUE), algorithm = "port",lower=c(0,0),upper=c(1,5))
summary(res)

plot(d$time, obs, xlab="Temps",ylab="H(t)", col="red")
lines(d$time, fitted(G),lwd=2, col="blue")


#Définir la fonction à minimiser pour utiliser l'algorithme génetique

Fitness<-function(x){
  par <- list(r=0.4,K=50000,m=0.03, c=30,omega=x[1], pi0=x[2],
              teta=80,m2=30,epsilon=0.1,epsilon2=0.1,
              f=4,pi2=0.9,a=x[3],sigma=0.9) # vecteurs de paramètres 
  V <- c(S=N-1, I=1,R=0,F=100 ,H=0)      # valeurs initiales du vecteur d'états
  c <- ode(V, 1:Tmax, deriv, parms=par)
  print(c[,"H"])
  s<- -mean((c[,"H"]-obs)^2)
  return(s)
}

# L'algorithme génétique 
G<-ga(type="real-value", fitness=Fitness, min=c(0,0,0), max=c(0.1,1,5),
  popSize=10, maxiter=10, pcrossover=0.5, pmutation=0, elitism=0)

plot(d$time, obs, xlab="Temps",ylab="H(t)", col="red")

#On dessine la courbe de l'évolution des infectés humais avec les paramètres estimés

N = 500                       # nombre d'élèves
Tmax = 120               # temps maximal
par <- list(r=0.4,K=50000,m=0.03, c=30,omega=0.05, pi0=0.81, teta=80,m2=30,epsilon=0.1,epsilon2=0.1,f=4,pi2=0.9,a=3.74,sigma=0.9) # vecteurs de paramètres 
V <- c(S=N-1, I=1,R=0,F=100 ,H=0)              # valeurs initiales du vecteur d'états
comparts2 <- ode(V, seq(1,Tmax,by=0.1), deriv, parms=par)
lines(d$time, comparts2[,"H"],lwd=2, col="blue")


#K_estim = coefficients(res)["K"]
# r_estim = coefficients(res)["r"]
# f_estim = coefficients(res)["f"]            
#Heq =   
#R0 = beta_estim/gamma_estim   # 1.24
# prédiction à Tmax+1 :
#H_pred = model_tofit(K_estim, pi0_estim, omega_estim, 1:(Tmax+1))[Tmax+1]
#H_pred = model_tofit(K_estim, pi0_estim, omega_estim, 1:(Tmax+1))[Tmax+1]

#Etude de sensibilité 
#Méthode de Sobol 
simulation <- function(paramMatrix){
  V <- c(S=N-1, I=1,R=0,F=100 ,H=0)   
  out <-NULL # j'initialise le vecteur renvoyé.
  for (i in 1:nrow(paramMatrix)){
    # simulation time 
    tps <- 1:Tmax
    # resolution de l'ODE pour la ième ligne de valeurs de paramètres :
    sol <- ode(V, tps, deriv, parms=c("r"=paramMatrix$r[i],"K"=paramMatrix$K[i],"m"=paramMatrix$m[i],"c"=paramMatrix$c[i],"omega"=paramMatrix$omega[i],"pi0"=paramMatrix$pi0[i],"teta"=paramMatrix$teta[i],"m2"=paramMatrix$m2[i],"epsilon"=paramMatrix$epsilon[i],"epsilon2"=paramMatrix$epsilon2[i],"f"=paramMatrix$f[i],"pi2"=paramMatrix$pi2[i],"a"=paramMatrix$a[i],"sigma"=paramMatrix$sigma[i]), method="lsoda")
    #print(sol[length(tps),"H"])
    out<-c(out, sol[Tmax,"H"])
    
  }
  return(out)
}


Nrep=100 # passer à plutot 10 000 une fois que c'est testé, pour avoir convergence.
# vérifier la convergence en lançant plusieurs fois la méthode Sobol et regardant si résultats similaires. 

# Matrice d'échantillonage :
X1<-data.frame(r=runif(Nrep, 0, 5), K=runif(Nrep, 0, 5),m=runif(Nrep, 0, 5),c=runif(Nrep, 0, 5),omega=runif(Nrep, 0, 5),pi0=runif(Nrep, 0, 5),teta=runif(Nrep, 0, 5),m2=runif(Nrep, 0, 5),epsilon=runif(Nrep, 0, 5),epsilon2=runif(Nrep, 0, 5),f=runif(Nrep, 0, 5),pi2=runif(Nrep, 0, 5),a=runif(Nrep, 0, 5),sigma=runif(Nrep, 0, 5))
# Matrice de ré-échantillonage :
X2<-data.frame(r=runif(Nrep, 0, 5), K=runif(Nrep, 0, 5),m=runif(Nrep, 0, 5),c=runif(Nrep, 0, 5),omega=runif(Nrep, 0, 5),pi0=runif(Nrep, 0, 5),teta=runif(Nrep, 0, 5),m2=runif(Nrep, 0, 5),epsilon=runif(Nrep, 0, 5),epsilon2=runif(Nrep, 0, 5),f=runif(Nrep, 0, 5),pi2=runif(Nrep, 0, 5),a=runif(Nrep, 0, 5),sigma=runif(Nrep, 0, 5))

# Faire un petit test de la fonction simulation, avant de passer à la suite.
simulation(X1)   



# indices d'ordre 1 et total
res2<-sobol2007(simulation, X1, X2, nboot=30)
# res2<-sobol2002(simulation, X1, X2, nboot=10)
plot(res2)
summary(res2)
print(res2)


## indices d'ordre 1 et 2 (facultatif)
# res<-sobol(simulation, X1, X2, order=2, nboot=10)
# plot(res)
# summary(res)
# print(res) 


# Pour méthode Morris : faire un print la 1ère fois permet de voir la
# matrice d'exploration de l'espace des paramètres (l'enlever ensuite).
simulationM <- function(paramMatrix){
  V <- c(S=N-1, I=1,R=0, F=100 ,H=0)
  out <-NULL
  
  for (i in 1:nrow(paramMatrix)){
    # parameters
    # ici param doit être au format liste (indifférent dans méthode précédente)
    #r=0.4,K=50000,m=0.03, c=30,omega=0.02, pi0=1, teta=80,m2=30,epsilon=0.1,epsilon2=0.1,f=4,pi2=0.9,a=4,sigma=0.9
    par <- list("r"=paramMatrix[i,1],"K"=paramMatrix[i,2],"m"=paramMatrix[i,3],"c"=paramMatrix[i,4],"omega"=paramMatrix[i,5],"pi0"=paramMatrix[i,6],"teta"=paramMatrix[i,7],"m2"=paramMatrix[i,8],"epsilon"=paramMatrix[i,9],"epsilon2"=paramMatrix[i,10],"f"=paramMatrix[i,11],"pi2"=paramMatrix[i,12],"a"=paramMatrix[i,13],"sigma"=paramMatrix[i,14]) 
    #print(unlist(par))
    tps <- seq(1, max(d$time),by=0.1)
    sol <- ode(V, tps, deriv, par, method="lsoda") 
    out<-c(out, sol[length(tps),"H"])
  }
  return(out)
}

resM <- morris(model = simulationM, factors = 14, r = 50, binf=0, bsup=1,
               design = list(type = "oat", levels = 10, grid.jump = 1))
plot(resM)
resM
