# Simulation mod�le de Code�o
# auteur AG. 01-2017.
rm(list=ls())

library(deSolve)
library(nls2)
library(sensitivity)

# ------------------ DONNEES
# Chemin du r�pertoire contenant les fichiers de donn�es :
path = "C:/Users/Augustin/Documents/Cours/Centrale/2A/Mod�lisation math�matique pour la biologie"
setwd(path)           
#Change le r�pertoire par d�faut

d1 <- read.table(paste(path, "/datab.txt", sep=""),
                 header=TRUE, sep="\t", na.strings="NA", dec=".")
d1
names(d1)=c("t","I","B1","D1","I2","B2","D2","C1","C2")

d2 <- read.table(paste(path, "/datab2.txt", sep=""),
                 header=TRUE, sep="\t", na.strings="NA", dec=".")
d2
names(d2)=c("t","I","B1","D1","I2","B2","D2","C1","C2")

d3 <- read.table(paste(path, "/datab3.txt", sep=""),
                 header=TRUE, sep="\t", na.strings="NA", dec=".")
d3
names(d3)=c("t","I","B1","D1","I2","B2","D2","C1","C2")

meteo <- read.table(paste(path, "/meteo.txt", sep=""),
                    header=TRUE, sep="\t", na.strings="NA", dec=".")


d=rbind(d1,d2,d3)

plot(d$t, d$I, xlab="temps", ylab="Nb infect�s")
plot(d$t, d$I2, xlab="temps", ylab="Nb infect�s")
plot(d$t, d$B1, xlab="temps", ylab="Concentration en bacteries")
plot(meteo$Day, meteo$Temp, xlab="temps", ylab="Concentration en bacteries",col="blue")
lines(meteo$Day, meteo$Temp, xlab="temps", ylab="Concentration en bacteries")

# ------------------ MODELE
# mod�lisation par compartiments : 
# V=vecteur d'�tats, p=vecteur de param�tres
deriv <- function(t,y,p)
{
  with(as.list(c(p, y)), {
    dS <- -a*(B/(K+B))*S*tan((1+m*sin((2*pi/g)*t)))  # mod�le de Code�o population constante et influence de l'environnement                                         # via taux de concentration B de toxig�ne
    dI <- a* (B/(K+B))*S*tan((1+m*sin((2*pi/g)*t)))-r*I
    dB <- -B*l+e*I                #l=nb-mb
    list(c(dS,dI,dB))
  })
}

# ------------------ SIMULATION

# test d'une simulation :
H = 1560                       # population totale du pays ou r�gion
Tmax = max(d$t)                 # temps maximal
par <- list( a=0.02, K=42000,r=0.17,l=0.04,e=2.9,m=0.5,g=40)    # vecteurs de param�tres, varient en fonction type �pid,pand,end�mie ... Cf feuille
V <- c(S=H-1, I=1, B=0)                            # valeurs initiales du vecteur d'�tats
comparts <- ode(V, 1:Tmax, deriv, parms=par)
plot(comparts[,"time"], comparts[,"I"],type="l",xlab="Temps",ylab="I(t)", lwd=4, col="red")  #temps en jour avec cette �quation

#test d'une simulation prolong�e en temps :
H = 1560                       # population totale du pays ou r�gion
Tmax = max(d$t)                 # temps maximal
resultat = testa
par <- list( a=0.07, K=52000,r=0.14,l=0.08,e=7,m=coefficients(resultat)["m"],g=coefficients(resultat)["g"])    # vecteurs de param�tres, varient en fonction type �pid,pand,end�mie ... Cf feuille
V <- c(S=H-1, I=1, B=0)                            # valeurs initiales du vecteur d'�tats
comparts <- ode(V, 1:(2*Tmax), deriv, parms=par)
plot(1:(2*Tmax), comparts[,"I"],type="l",xlab="Temps",ylab="I(t)", lwd=4, col="red")  #temps en jour avec cette �quation

  # Fonction pour l'estimation des param�tres
  # Retourne un vecteur I simul�, qui doit �tre de m�me taille que celui des donn�es.
model_tofit <- function(a,K,r,l,e,m,g, times){
  V <- c(S=H-1, I=1,B=0)
  Fsim <- ode(V, times, deriv, parms=list( "a"=a,"K"=K,"r"=r,"l"=l,"e"=e,"m"=m,"g"=g))
  #print(Fsim)
  return(Fsim[,"I"])  # note : il faut que le vecteur times d�marre � 1, sinon le temps initial est incorrect
                      # sinon mettre c(1, times) au lieu de times dans ode et enlever le 1er terme du vecteur renvoy�.
}
  # estimation des param�tres p :
obs <- d$I
#res<-nls2(obs~model_tofit(a,K,r,l,e,m, T), start=list(m=0.0001), data=list(a=1,r=0.02,K=100000, l=0.04,e=10,T=d$t), trace=TRUE, nls.control(warnOnly=TRUE), algorithm = "port")  # pour garder des param�trse constants : data=list(beta=1),
#res2<-nls2(obs~model_tofit(a,K,r,l,e,m, T), start=data.frame( a=c(0, 1)), data=list(T=d$t, K=100000, e=10, r=0.02,l=0.03, m=0.0001), trace=TRUE, nls.control(warnOnly=TRUE, maxiter=10), algorithm = "port", all=FALSE)
#res3<-nls2(obs~model_tofit(a,K,r,l,e,m, T), start=list(a=1,r=0.2,l=0.33,e=10,m=1), data=list(K=100000,T=d$t), trace=TRUE, nls.control(warnOnly=TRUE), algorithm = "port")
#res4<-nls2(obs~model_tofit(a,K,r,l,e,m, T), start=list(a=1, K=100000,r=0.02,l=0.04,e=10, m=0.0001), data=list(T=d$t), trace=TRUE, nls.control(warnOnly=TRUE), algorithm = "port",lower=c(0.00001,100000,0.07,0.02,0.01,0), upper=c(1,1000000,0.4,0.4,10,1))
#result<-nls2(obs~model_tofit(a,K,r,l,e,m, T), start=data.frame( a=c(0.00001, 1), K=c(100000,1000000),r=c(0.07,0.4),l=c(0.02,0.4),e=c(0.01,10),m=c(0.00001,100),g=c(1,100)), data=list(T=d$t), trace=TRUE, nls.control(warnOnly=TRUE, maxiter=10), algorithm = "port")
testo<-nls2(obs~model_tofit(a,K,r,l,e,m,g, T), start=data.frame(m=c(0.01,1),g=c(50,80)), data=list(a=0.07,K=52000,r=0.14,l=0.08,e=7,T=d$t), trace=TRUE, nls.control(warnOnly=TRUE, maxiter=10), algorithm = "port")
testa<-nls2(obs~model_tofit(a,K,r,l,e,m,g, T), start=data.frame( a=c(0.1, 0.3), K=c(10000,100000),r=c(0.05,0.4),l=c(0.01,0.2),e=c(0.1,10),m=c(0.0001,10),g=c(40,80)), data=list(T=d$t), trace=TRUE, nls.control(warnOnly=TRUE, maxiter=10), algorithm = "port")
testi<--nls2(obs~model_tofit(a,K,r,l,e,m,g, T), start=data.frame( a=c(0.1, 0.3), K=c(10000,100000),r=c(0.05,0.4),l=c(0.01,0.2),e=c(0.1,10),m=c(0.0001,10)), data=list(g=40,T=d$t), trace=TRUE, nls.control(warnOnly=TRUE, maxiter=10), algorithm = "port")

  # fixer les param�tres donn�es par l'enonc� et faire varier ceux inconnus
  # data sert pour garder des param�trse constants 
  # warn.only=TRUE permet de stocker le dernier r�sultat obtenu m�me si l'algo �choue
  # l'utilisation de l'algo "port" permet de fixer des bornes � l'espace de recherche des param�tres (ici on veut les garder positifs par exemple)
  # nls2 permet de relancer plusieurs fois l'optimisation en partant de diff�rents points initiaux
  # (note : algos brute-force et grid-search ou random-search font seulement des �chantillonages de la fonction objectif)

"summary(res)
summary(res2)
summary(res3)
summary(res4)
summary(result)"
summary(testo)
summary(testa)
summary(testi)

  # Parameters:
  #   Estimate Std. Error t value Pr(>|t|)    
  # gamma   4.9028     0.4424   11.08 0.000104 ***
  #   beta    6.0887     0.5113   11.91 7.36e-05 ***
  #   ---
  #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
  # 
  # Residual standard error: 0.4125 on 5 degrees of freedom

plot(d$t, obs, xlab="Temps",ylab="I(t)", col="red") 
#lines(d$t, fitted(res),lwd=2, col="blue")
#lines(d$t, fitted(res2),lwd=2, col="green")
#lines(d$t, fitted(res3),lwd=2, col="pink")
#lines(d$t, fitted(res4),lwd=2, col="purple")
#lines(d$t, fitted(result),lwd=2, col="yellow")
lines(d$t, fitted(testo),lwd=2, col="yellow")
lines(d$t, fitted(testa),lwd=2, col="black")
lines(d$t, fitted(testi),lwd=2, col="pink")

# # facultatif : pr�diction et confiance (bugs sur certaines versions...)
# confint(res)  # intervalles de confiance sur les estim�es des param�tres
# conf<-predict(as.lm(res), interval="confidence")
# lines(d$t, conf[,"lwr"], col="green")
# lines(d$t, conf[,"upr"], col="green")
# conf<-predict(as.lm(res), interval="prediction")
# lines(d$t, conf[,"lwr"], col="dark green")
# lines(d$t, conf[,"upr"], col="dark green")

# ------------------ UTILISATION DES VALEURS ESTIMEES
res3=testo
a_est = coefficients(res3)["a"]
K_est = coefficients(res3)["K"]
r_est = coefficients(res3)["r"]
l_est = coefficients(res3)["l"]
e_est = coefficients(res3)["e"]
m_est = coefficients(res3)["m"]
g_est = coefficients(res3)["g"]
Ieq = N*(1-gamma_estim/beta_estim)  #9.35
R0 = (H*e_est*a_est)/(l_est*r_est*K_est)   # 1.24
  # pr�diction � Tmax+1 :
I_pred = model_tofit(a_est,K_est,r_est,l_est,e_est,m_est,g_est, 1:(Tmax+24))[Tmax+24]
I_pred  


# ------------------ 2�) que se passe-t-il si on n'utilise que les premi�res donn�es ?
# Par exemple t<= Tmax/2
Thalf = round(Tmax/2)
obs2 <- d$I[d$t<=Thalf]
timeline <- d$t[d$t<=Thalf]
res2<-nls(obs2~model_tofit(beta, gamma, T), start=list(gamma=0.2, beta=1), data=list(T=timeline), trace=TRUE)  # pour garder des param�trse constants : data=list(beta=1), 
plot(timeline, obs2, xlab="Temps",ylab="I(t)", col="red") 
lines(timeline, fitted(res2),lwd=2, col="blue")
  # gamma=2.498  beta= 3.470 : pas les m�mes valeurs qu'avec le jeu complet.      
R0 = coefficients(res2)["beta"]/coefficients(res2)["gamma"] 

# ------------------ 3�) que se passe-t-il si on utilise les donn�es estim�es � partir de l'�chantillonage ?
Nech=23
plot(d$t, d$Iech/Nech, ylab = "I/N", col="red", ylim=range(0, 1.5*d$Iech/Nech,na.rm=TRUE ))
points(d$t, d$I/N, col="black")
  # intervalle de confiance � 95% :
  # conditions : f=N_i/Nech n'est pas proche de 0 ou de 1 ; N_i >= 5 ; Nech-Ni>=5
  # (conditions pas toujours v�rifi�es !!)
q=qnorm(1-0.05/2)
delta = 2*q*sqrt(d$Iech/Nech*(1-d$Iech/Nech)/Nech)
Ilow = d$Iech/Nech - delta/2 
Iup  = d$Iech/Nech + delta/2 
lines(d$t, Ilow, col="red")
lines(d$t, Iup, col="red")
# tr�s larges intervalles de confiance, comme attendu.
# puis refaire calculs avec Iech/Nech *N comme donn�es au lieu de I.
# Ne garder que les valeurs non NA :
obs3 <- d$Iech*N/Nech
timeline3 <- d$t[which(!is.na(obs3))]
obs3 <- obs3[!is.na(obs3)]
res3<-nls(obs3~model_tofit(beta, gamma, T), start=list(gamma=4, beta=6), data=list(T=timeline3), trace=TRUE, nls.control(warnOnly=TRUE))   
summary(res3)
plot(timeline3, obs3, xlab="Temps",ylab="I(t)", col="red") 
lines(timeline3, fitted(res3),lwd=2, col="blue")
# En enlevant la derni�re mesure, peut converger mais erreurs standard �normes !
Ieq = N*(1-coefficients(res3)["gamma"]/coefficients(res3)["beta"])  
R0 = coefficients(res3)["beta"]/coefficients(res3)["gamma"]  

# interpr�tation biologique : 
#  en r�gime permanent : la proba de gu�rir correspond � celle de transmettre le virus � un voisin
# mais le probl�me est qu'il y a la phase d'�mergence.
# cette phase fait aussi que la dur�e moyenne d'infection est biais�e.
# La premi�re personne infect�e le reste pendant une dur�e tr�s longue, 
# qui diminue pour les personnes suivantes.

# ------------------  4�) Analyse de sensibilit� globale du mod�le
# fonction qui prend en entr�e un "data-frame" de param�tres 
# (les X[i,p] o� 1<=i<=n o� n nombre de param�tres ; 1<=p<=Nrep o� Nrep nombre de r�p�titions)
# et retourne un vecteur Yp de sorties de taille Nrep.
# (une sortie par ligne du data-frame)
simulation <- function(paramMatrix){
  V <- c(S=H-65357, I=65357,B=0)   
  out <-NULL # j'initialise le vecteur renvoy�
  # on parcourt chacune des lignes du data-frame : une ligne = 1 valeur d'un jeux de param�tres.
  for (i in 1:nrow(paramMatrix)){
    # simulation time 
    tps <- seq(1, max(d$t),by=0.5)
    # resolution de l'ODE pour la i�me ligne de valeurs de param�tres :
    sol <- ode(V, tps, deriv, parms=c("a"=paramMatrix$a[i],"K"=paramMatrix$K[i],"r"=paramMatrix$r[i],"l"=paramMatrix$l[i],"e"=paramMatrix$e[i]), method="lsoda")
    # Remarque : on peut aussi utiliser la m�thode rk4 pour gagner en pr�cision, si besoin. Voir l'aide de 'ode'.
    # ici j'ai choisi comme variable d'int�r�t le nombre d'infect�s au temps final
    # mais on pourrait imaginer d'autres variables comme par exemple la somme des carr�s
    # utilis�es pour model_tofit pour l'estimation :
    out<-c(out, sol[length(tps),"I"])
  }
  return(out)
}

Nrep=10000 # passer � plutot 10 000 une fois que c'est test�, pour avoir convergence.
         # v�rifier la convergence en lan�ant plusieurs fois la m�thode Sobol et regardant si r�sultats similaires. 

  # Matrice d'�chantillonage :
X1<-data.frame(a=runif(Nrep,0.00001,1),K=runif(Nrep,100000,1000000),r=runif(Nrep,0.07,0.4),l=runif(Nrep,0.02,0.4),e=runif(Nrep,0.02,10))
  # Matrice de r�-�chantillonage :
X2<-data.frame(a=runif(Nrep,0.00001,1),K=runif(Nrep,100000,1000000),r=runif(Nrep,0.07,0.4),l=runif(Nrep,0.02,0.4),e=runif(Nrep,0.02,10))
# remarque: on peut pour certains param�tres mettre une valeur constante, 
# ce qui revient � consid�rer qu'ils sont connus avec pr�cision et 
# donc leur incertitude n'est pas prise en compte dans le calcul des indices.

# Faire un petit test de la fonction simulation, avant de passer � la suite.
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

# Explications :
# --------------
# La fonction d�finissant le mod�le pour Sobol prend en argument une matrice de taille Nrep x p dont chaque ligne repr�sente une r�alisation 
# d'un jeu de p param�tres, tir�e selon les lois d�finies pour caract�riser l'incertitude de chacun d'entre eux. A l'int�rieur de cette fonction,
# on va lancer autant de simulations qu'on donne de jeux de param�tres diff�rents, c'est-�-dire Nrep, le nombre de lignes de la matrice (nrow()).
# Pour chaque simulation, il faut d�finir une sortie scalaire (ie ce doit �tre un r�el : si on en veut plusieurs, il faut lancer autant de fois l'analyse). 
# Dans mon exemple, j'ai choisi de prendre la valeur du nombre d'infect�s au dernier temps de ma simulation (Tmax). A vous de choisir une sortie d'int�r�t
# pour votre mod�le. La fonction doit renvoyer un vecteur dont chaque composante correspond � la sortie d'une simulation : ce vecteur est donc de taille Nrep.
# Ensuite la fonction Sobol2002/7 prend en argument le nom de cette fonction, et deux matrices de jeux de param�tres X1 et X2. Elle comprend aussi un argument de
# nombre de r�p�tition de "bootstrap" qui sert � estimer l'intervalle de confiance sur la valeur de l'indice.

# Pour m�thode Morris : faire un print la 1�re fois permet de voir la
# matrice d'exploration de l'espace des param�tres (l'enlever ensuite).
simulationM <- function(paramMatrix){
  V <- c(S=N-1, I=1) 
  out <-NULL
  
  for (i in 1:nrow(paramMatrix)){
    # parameters
    # ici param doit �tre au format liste (indiff�rent dans m�thode pr�c�dente)
    par <- list("beta"=paramMatrix[i,1],"gamma"=paramMatrix[i,2]) 
    #print(unlist(par))
    tps <- seq(1, max(d$t),by=0.1)
    sol <- ode(V, tps, deriv, par, method="lsoda") 
    out<-c(out, sol[length(tps),"I"])
  }
  return(out)
}
resM <- morris(model = simulationM, factors = 2, r = 50, binf=0, bsup=1,
               design = list(type = "oat", levels = 10, grid.jump = 1))
plot(resM)
resM


