library(deSolve)
library(nls2)
library(sensitivity)



# ------------------ DONNEES
# Chemin du répertoire contenant les fichiers de données :
path = "C:/Users/Imane_Benhayoun/Desktop/Docs_cours_centrale/2ème_année/Maths_pour_la_bio"
setwd(path)                      #Change le répertoire par défaut

d <- read.table(paste(path, "/data_class_2018.txt", sep=""),
                header=TRUE, sep="\t", na.strings="NA", dec=".")
d
plot(d$t, d$I, xlab="temps", ylab="Nb infectés")


# ------------------ MODELE
# modélisation par compartiments : 
# V=vecteur d'états, p=vecteur de paramètres

deriv <- function(t, V, parametre)
{
  # Janvier 2015
  #Fonction donnant l'evolution du modele avec deux SEIR
  
  with (as.list(c(parametre, V)), {
    lambda_A = a_A *( beta_A * I_A + beta_B * I_B )
    lambda_B = a_B *( beta_A * I_A + beta_B * I_B )
    
    dS_A <- - lambda_A * S_A
    dE_A <- lambda_A * S_A - delta_SEIR * E_A
    dI_A <- delta_SEIR *E_A - gamma *I_A
    dR_A <- gamma * I_A
    dS_B <- - lambda_B * S_B
    dE_B <- lambda_B * S_B - delta_SEIR * E_B
    dI_B <- delta_SEIR *E_B - gamma *I_B
    dR_B <- gamma * I_B
    list(c(dS_A,dE_A,dI_A,dR_A,dS_B,dE_B,dI_B,dR_B))
  })
}

# ------------------ SIMULATION
a_A = 6
beta_A = 8.86
a_B = 3
beta_B = 4.43
delta_SEIR = 100
gamma = 74
Tmax = 10
N= 10000  #Nombre total d'individus
parametre <- list(a_A , beta_A , a_B , beta_B , delta_SEIR , gamma)
V <- c(S_A=0.82, E_A=0, I_A=0.06, R_A=0, S_B=0.11, E_B=0, I_B=0.01, R_B=0 )
TT = seq(0,1,length=100)
comparts <- ode(V, TT, deriv, parms=parametre)

bruit <- 0.005*rnorm(Tmax,mean=0,sd=0.2)


d <- data.frame(time=TT, S_A=(comparts[,"S_A"]) , S_B=(comparts[,"S_B"]), E_A=(comparts[,"E_A"]) , E_B=(comparts[,"E_B"]), I_A=(comparts[,"I_A"]) , I_B=(comparts[,"I_B"]), R_A=(comparts[,"R_A"]) , R_B=(comparts[,"R_B"]))

#par(mfrow=c(2,2), mar=c(4,5,0.5,1))
plot(comparts[,"time"], comparts[,"S_A"],type="l", xlab="Temps",ylab="S_A(t)", lwd=4, col="red")
plot(comparts[,"time"], comparts[,"E_A"],type="l", xlab="Temps",ylab="E_A(t)", lwd=4, col="blue")
plot(comparts[,"time"], comparts[,"I_A"],type="l", xlab="Temps",ylab="I_A(t)", lwd=4, col="green")
plot(comparts[,"time"], comparts[,"R_A"],type="l", xlab="Temps",ylab="R_A(t)", lwd=4, col="black")

plot(comparts[,"time"], comparts[,"S_B"],type="l", xlab="Temps",ylab="S_B(t)", lwd=4, col="red")
plot(comparts[,"time"], comparts[,"E_B"],type="l", xlab="Temps",ylab="E_B(t)", lwd=4, col="blue")
plot(comparts[,"time"], comparts[,"I_B"],type="l", xlab="Temps",ylab="I_B(t)", lwd=4, col="green")
plot(comparts[,"time"], comparts[,"R_B"],type="l", xlab="Temps",ylab="R_B(t)", lwd=4, col="black")

par(mfrow=c(2,2), mar=c(4,5,0.5,1))
plot(comparts[,"time"], comparts[,"S_A"]+bruit,type="l", xlab="Temps",ylab="S_A(t)_bruit", lwd=2, col="red")
plot(comparts[,"time"], comparts[,"E_A"]+bruit,type="l", xlab="Temps",ylab="E_A(t)_bruit", lwd=2, col="blue")
plot(comparts[,"time"], comparts[,"I_A"]+bruit,type="l", xlab="Temps",ylab="I_A(t)_bruit", lwd=2, col="green")
plot(comparts[,"time"], comparts[,"R_A"]+bruit,type="l", xlab="Temps",ylab="R_A(t)_bruit", lwd=2, col="black")

plot(comparts[,"time"], comparts[,"S_B"]+bruit,type="l", xlab="Temps",ylab="S_B(t)_bruit", lwd=2, col="red")
plot(comparts[,"time"], comparts[,"E_B"]+bruit,type="l", xlab="Temps",ylab="E_B(t)_bruit", lwd=2, col="blue")
plot(comparts[,"time"], comparts[,"I_B"]+bruit,type="l", xlab="Temps",ylab="I_B(t)_bruit", lwd=2, col="green")
plot(comparts[,"time"], comparts[,"R_B"]+bruit,type="l", xlab="Temps",ylab="R_B(t)_bruit", lwd=2, col="black")


# ------------------ ESTIMATION
# Fonction pour l'estimation des paramètres
# Retourne des vecteurs simulés , qui doit être de même taille que celui des données.

model_tofit <- function(a_A , a_B, beta_A , beta_B, delta_SEIR , gamma, times){
  V <- c(S_A=0.82, E_A=0, I_A=0.06, R_A=0, S_B=0.11, E_B=0, I_B=0.01, R_B=0 )
  Fsim <- ode(V, times, deriv, parms=list('a_A'=a_A , 'a_B'=a_B, 'beta_A'=beta_A , 'beta_B'=beta_B,  'delta_SEIR'=delta_SEIR , 'gamma'=gamma))
  return(c(Fsim[,"S_A"],Fsim[,"E_A"], Fsim[,"I_A"], Fsim[,"R_A"], Fsim[,"S_B"], Fsim[,"E_B"], Fsim[,"I_B"], Fsim[,"R_B"]) )
}

model_tofit_bruit <- function(a_A , a_B, beta_A , beta_B, delta_SEIR , gamma, times){
  V <- c(S_A=0.82, E_A=0, I_A=0.06, R_A=0, S_B=0.11, E_B=0, I_B=0.01, R_B=0 )
  Fsim <- ode(V, times, deriv, parms=list('a_A'=a_A , 'a_B'=a_B, 'beta_A'=beta_A , 'beta_B'=beta_B,  'delta_SEIR'=delta_SEIR , 'gamma'=gamma))
  return(c(Fsim[,"S_A"]+bruit,Fsim[,"E_A"]+bruit, Fsim[,"I_A"]+bruit, Fsim[,"R_A"]+bruit, Fsim[,"S_B"]+bruit, Fsim[,"E_B"]+bruit, Fsim[,"I_B"]+bruit, Fsim[,"R_B"]+bruit) )
}


obs <- model_tofit_bruit(a_A , a_B, beta_A , beta_B, delta_SEIR , gamma, TT)

res<- nls2(obs~model_tofit(a_A , a_B, beta_A , beta_B, delta_SEIR , gamma, TT), start=list(a_A=25, a_B=50, beta_A=10, beta_B=5, delta_SEIR=80, gamma=4), data=list(T=TT), trace=TRUE, nls.control(warnOnly=TRUE), algorithm = "port", lower=c(0, 0))  

summary(res)

L <- list("Estimation S_A(t)","Estimation E_A(t)","Estimation I_A(t)","Estimation R_A(t)","Estimation S_B(t)","Estimation E_B(t)","Estimation I_B(t)","Estimation R_B(t)")

#par(mfrow=c(2,2), mar=c(4,5,0.5,1) )

for (i in 1:8){
  a = length(TT)*(i-1)+1
  b = length(TT)*i
  plot(TT, fitted(res)[a:b],xlab="Temps",ylab=L[i], type="l", lwd=2)}

# ------------------ PREDICTION

a_A_estim = coefficients(res)["a_A"]
a_B_estim = coefficients(res)["a_B"]
beta_A_estim = coefficients(res)["beta_A"]
beta_B_estim = coefficients(res)["beta_B"]
delta_SEIR_estim = coefficients(res)["delta_SEIR"]
gamma_estim = coefficients(res)["gamma"]
# Ieq = N-gamma_estim/beta_estim 
# R0 = beta_estim*N/gamma_estim

# prédiction à Tmax+1 :
I_pred = model_tofit(a_A_estim, a_B_estim, beta_A_estim, beta_B_estim, delta_SEIR_estim, gamma_estim , 1:(length(TT)+1))[length(TT)+1]

par(mfrow=c(1,1), mar=c(4,5,1.5,1.5) )

# ------------------ Analyse de sensibilité globale du modèle
# Pour méthode Morris : faire un print la 1ère fois permet de voir la
# matrice d'exploration de l'espace des paramètres (l'enlever ensuite).
simulationM <- function(paramMatrix){
  V <- c(S_A=0.82, E_A=0, I_A=0.06, R_A=0, S_B=0.11, E_B=0, I_B=0.01, R_B=0 ) 
  out <-NULL
  for (i in 1:nrow(paramMatrix)){
    # parameters
    # ici param doit être au format liste (indifférent dans méthode précédente)
    par <- list("a_A" = paramMatrix[i,1], "a_B" = paramMatrix[i,2], "beta_A" = paramMatrix[i,3], "beta_B" = paramMatrix[i,4], "delta_SEIR"=paramMatrix[i,5],"gamma"=paramMatrix[i,6]) 
    #print(unlist(par))
    tps <- TT
    sol <- ode(V, tps, deriv, par, method="lsoda") 
    out<-c(out, sol[length(tps),"I_A"])
  }
  return(out)
}
resM <- morris(model = simulationM, factors = 6, r = 50, binf=0, bsup=1,
               design = list(type = "oat", levels = 10, grid.jump = 1))
plot(resM)
resM



#SOBOL ne marche pas je ne sais pas pourquoi ..

# fonction qui prend en entrée un "data-frame" de paramètres 
# (les X[i,p] où 1<=i<=n où n nombre de paramètres ; 1<=p<=Nrep où Nrep nombre de répétitions)
# et retourne un vecteur Yp de sorties de taille Nrep.
# (une sortie par ligne du data-frame)
#simulation <- function(paramMatrix){
#  V <- c(S_A=0.82, E_A=0, I_A=0.06, R_A=0, S_B=0.11, E_B=0, I_B=0.01, R_B=0 )  
#  out <-NULL # j'initialise le vecteur renvoyé
#  # on parcourt chacune des lignes du data-frame : une ligne = 1 valeur d'un jeux de paramètres.
#  for (i in 1:nrow(paramMatrix)){
#    # simulation time 
#    tps <- TT 
#    # resolution de l'ODE pour la ième ligne de valeurs de paramètres :
#    sol <- ode(V, tps, deriv, parms=c("a_A" = paramMatrix$a_A[i], "a_B" = paramMatrix$a_B[i], "beta_A" = paramMatrix$beta_A[i], "beta_B" = paramMatrix$beta_B[i], "delta_SEIR"=paramMatrix$delta_SEIR[i],"gamma"=paramMatrix$gamma[i], method="lsoda"))
#    # Remarque : on peut aussi utiliser la méthode rk4 pour gagner en précision, si besoin. Voir l'aide de 'ode'.
#    # ici j'ai choisi comme variable d'intérêt le nombre d'infectés au temps final
#    # mais on pourrait imaginer d'autres variables comme par exemple la somme des carrés
#    # utilisées pour model_tofit pour l'estimation :
#    out<-c(out, sol[length(tps),"I_A"])
#  }
#  return(out)
#}

#Nrep=100 # passer à plutot 10 000 une fois que c'est testé, pour avoir convergence.
# vérifier la convergence en lançant plusieurs fois la méthode Sobol et regardant si résultats similaires. 

# Matrice d'échantillonage :
#X1<-data.frame(a_A=runif(Nrep, 0, 10), a_B=runif(Nrep, 0, 10), beta_A=runif(Nrep, 0, 10), beta_B=runif(Nrep, 0,10), delta_SEIR=runif(Nrep, 0, 100), gamma=runif(Nrep, 0, 100))
# Matrice de ré-échantillonage :
#X2<-data.frame(a_A=runif(Nrep, 0, 10), a_B=runif(Nrep, 0, 10), beta_A=runif(Nrep, 0, 10), beta_B=runif(Nrep, 0,10), delta_SEIR=runif(Nrep, 0, 100), gamma=runif(Nrep, 0, 100))
# remarque: on peut pour certains paramètres mettre une valeur constante, 
# ce qui revient à considérer qu'ils sont connus avec précision et 
# donc leur incertitude n'est pas prise en compte dans le calcul des indices.

# Faire un petit test de la fonction simulation, avant de passer à la suite.
#simulation(X1)

# indices d'ordre 1 et total
#res2<-sobol2007(simulation, X1, X2, nboot=30)
#res2<-sobol2002(simulation, X1, X2, nboot=10)
#plot(res2)
#summary(res2)
#print(res2)
