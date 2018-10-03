# Simulation de l'exercice propagation d'épidémie dans la classe.
# auteur VL. 
rm(list=ls())

library(deSolve)
library(nls2)    # pas obligatoire (on peut utiliser nls au lieu de nls2)
library(sensitivity)



# ------------------ DONNEES
# Chemin du répertoire contenant les fichiers de données :
path = "C:/Users/ElAamrani/Desktop"
setwd(path)                      #Change le répertoire par défaut

d <- read.table(paste(path, "/data_class_2018.txt", sep=""),
                  header=TRUE, sep="\t", na.strings="NA", dec=".")
d
plot(d$t, d$I, xlab="temps", ylab="Nb infectés")




# ------------------ MODELE
  # modélisation par compartiments : 
  # V=vecteur d'états, p=vecteur de paramètres
# deriv <- function(t,y,p)
# {
#   with(as.list(c(p, y)), {
#     dS <- (-beta*S*I+gamma*I)*t
#     dI <- (beta*S*I-gamma*I)*t
#     list(c(dS,dI))
#   })
# }

# # ODE SIR-SI function
# SIRSI <- function(times, state, parameters) {
#   
#   with(as.list(c(state, parameters)), {
#     
#     Nh <- Sh + Ih + Rh
#     
#     dSh <-  (-(betah * bitting)/Nh) * (Sh * Iv)
#     dIh <-  ((betah * bitting)/Nh) * (Sh * Iv) - (gamma * Ih)
#     dRh <-  (gamma * Ih)
#     
#     Nv <- Sv + Iv
#     
#     dSv <-  (-((betav * bitting)/Nh) * (Sv * Ih))
#     dIv <-  (((betav * bitting)/Nh) * (Sv * Ih))
#     
#     return(list(c(dSh, dIh, dRh, dSv, dIv), Nh = Nh, Nv = Nv))
#   })
# }


SIRSI <- function(times, state, parameters) {

  with(as.list(c(state, parameters)), {

    Nh <- Sh + Ih + Rh

    dSh <-  muh * Nh - (muh + p + Cvh * Iv/Nh) * Sh
    dIh <-  (Cvh * Iv/Nh) * Sh - (gammah + muh) * Ih
    dRh <-  p * Sh + gammah * Ih - muh * Rh

    Nv <- Sv + Iv

    dSv <-  muv * Nv - (muv + Chv * Ih/Nh) * Sv
    dIv <-  (Chv * Ih/Nh) * Sv - muv * Iv

    return(list(c(dSh, dIh, dRh, dSv, dIv), Nh = Nh, Nv = Nv))
  })
}

  
# ------------------ SIMULATION
# test d'une simulation :
N = 10001                          # nombre de susceptible
Tmax = 50                  # temps maximal
# par <- list(betah=0.4, betav=0.4, bitting=1, gamma=0.167) # vecteurs de paramètres 
par <- list(muh=1/25000, muv=1/4, p=0, Cvh=0.75, Chv=0.375, gammah=1/3-1/25000)
par1 <- list(muh=1/25000, muv=1/4, p=0.1, Cvh=0.75, Chv=0.375, gammah=1/3-1/25000)
par2 <- list(muh=1/25000, muv=1/4, p=0.2, Cvh=0.75, Chv=0.375, gammah=1/3-1/25000)
par3 <- list(muh=1/25000, muv=1/4, p=0.3, Cvh=0.75, Chv=0.375, gammah=1/3-1/25000)
state <- c(Sh=10000, Ih=10, Rh=0, Sv=50000, Iv=10)              # valeurs initiales du vecteur d'états
comparts <- ode(state, 1:Tmax, SIRSI, parms=par)
comparts1 <- ode(state, 1:Tmax, SIRSI, parms=par1)
comparts2 <- ode(state, 1:Tmax, SIRSI, parms=par2)
comparts3 <- ode(state, 1:Tmax, SIRSI, parms=par3)

par(mfrow=c(1,2))
plot(comparts[,"time"], comparts[,"Ih"],type="l",xlab="Temps",ylab="Ih(t)", lwd=4, col="red", sub = "Evolution des infectés dans la population humaine", pch = 1)
lines(comparts[,"time"], comparts1[,"Ih"],type="l",xlab="Temps",ylab="Ih(t)", lwd=4, col="pink")
lines(comparts[,"time"], comparts2[,"Ih"],type="l",xlab="Temps",ylab="Ih(t)", lwd=4, col="blue")
lines(comparts[,"time"], comparts3[,"Ih"],type="l",xlab="Temps",ylab="Ih(t)", lwd=4, col="green")
legend(25,3000, legend=c("p=0","p=0.1","p=0.2","p=0.3"),lty=c(1,1,1,1), cex=0.8,
    lwd=c(2.5,2.5,2.5,2.5),col=c("red","pink","blue","green")) # gives the legend lines the correct color and width

plot(comparts[,"time"], comparts[,"Iv"],type="l",xlab="Temps",ylab="Iv(t)", lwd=4, col="red", sub =  "Evolution des infectés dans la population des moustiques", pch = 2)
lines(comparts[,"time"], comparts1[,"Iv"],type="l",xlab="Temps",ylab="Ih(t)", lwd=4, col="pink")
lines(comparts[,"time"], comparts2[,"Iv"],type="l",xlab="Temps",ylab="Ih(t)", lwd=4, col="blue")
lines(comparts[,"time"], comparts3[,"Iv"],type="l",xlab="Temps",ylab="Ih(t)", lwd=4, col="green")
legend(25,10400, legend=c("p=0","p=0.1","p=0.2","p=0.3"),lty=c(1,1,1,1), cex=0.8,
       lwd=c(2.5,2.5,2.5,2.5),col=c("red","pink","blue","green"))




# ------------------ ESTIMATION
# Fonction pour l'estimation des paramètres
# Retourne un vecteur I simulé, qui doit être de même taille que celui des données.
model_tofit <- function(muh, muv, p, Cvh, Chv, gammah, times){
  state <- c(Sh=10000, Ih=1, Rh=0, Sv=50000, Iv=1) 
  Fsim <- ode(state, times, SIRSI, parms=list("muh"=muh, "muv"=muv, "p"=p, 
                                        "Cvh"=Cvh, "Chv"=Chv, "gammah"=gammah))
  return(Fsim[,"Ih"])  
}


model_tofit2 <- function(muh, muv, p, Cvh, Chv, gammah, times){
  state <- c(Sh=10000, Ih=1, Rh=0, Sv=50000, Iv=1) 
  Fsim <- ode(state, times, SIRSI, parms=list("muh"=muh, "muv"=muv, "p"=p, 
                                        "Cvh"=Cvh, "Chv"=Chv, "gammah"=gammah))
  return(c(Fsim[,"Ih"], Fsim[,"Iv"]))  
      # note : il faut que le vecteur times démarre à 1, sinon le temps initial est incorrect
      # sinon mettre c(1, times) au lieu de times dans ode et enlever le 1er terme du vecteur renvoyé.
}


#------------------------------------------------------------------------------
# Estimation basique
T2 = Tmax
obs <- comparts[,"Ih"][1:T2] + rnorm(T2, mean=0, sd=2)
res<-nls2(obs~model_tofit(muh, muv, p, Cvh, Chv, gammah, T), 
          start=list(muh = 0.00003, muv = 0.2, Cvh=0.8, Chv=0.4, gammah=0.3), 
          d=list(T=1:T2, p=0),
          trace=TRUE, 
          nls.control(warnOnly=TRUE), 
          algorithm = "port",
          lower = c(0.000027,0.1,0.6,0.2,0.2), upper = c(0.00014,0.4,0.8,0.5,0.5))
# d=list() sert pour garder des paramètrse constants 
# warn.only=TRUE permet de stocker le dernier résultat obtenu même si l'algo échoue
# l'utilisation de l'algo "port" permet de fixer des bornes à l'espace de recherche des paramètres (ici on veut les garder positifs par exemple)
# nls2 permet de relancer plusieurs fois l'optimisation en partant de différents points initiaux
# (note : algos brute-force et grid-search ou random-search font seulement des échantillonages de la fonction objectif)


summary(res)
par(mfrow=c(1,2))



N = 10001                          # nombre de susceptible
Tmax = 50                  # temps maximal
# par <- list(betah=0.4, betav=0.4, bitting=1, gamma=0.167) # vecteurs de paramètres 
par <- list(muh=coefficients(res)["muh"], muv=coefficients(res)["muv"],
            p=0, Cvh=coefficients(res)["Cvh"], Chv=coefficients(res)["Chv"], 
            gammah=coefficients(res)["gammah"])
state <- c(Sh=10000, Ih=1, Rh=0, Sv=50000, Iv=1)              # valeurs initiales du vecteur d'états
comparts2 <- ode(state, 1:Tmax, SIRSI, parms=par)


plot(comparts[,"time"], comparts2[,"Ih"],type="o",xlab="Temps",ylab="Ih(t)", col="blue",  sub="Prédiction dans la population humaine", pch = 1)
lines(comparts[,"time"], comparts[,"Ih"],type="l",xlab="Temps",ylab="Ih(t)", col="red")

plot(comparts[,"time"], comparts2[,"Iv"],type="o",xlab="Temps",ylab="Iv(t)", col="blue",  sub="Prédiction dans la population des moustiques", pch = 2)
lines(comparts[,"time"], comparts[,"Iv"],type="l",xlab="Temps",ylab="Iv(t)", col="red")










#-------------------------------------------------------------------------------
# estimation des paramètres p :
obs <- comparts[,"Ih"][1:T2] + rnorm(T2, mean=0, sd=2)
obs2 <- comparts[,"Iv"][1:T2] + rnorm(T2, mean=0, sd=2)
res2<-nls2(c(obs, obs2)~model_tofit2(muh, muv, p, Cvh, Chv, gammah, T), 
          start=list(gammah=0.3, Cvh=0.8, muv = 0.2, Chv=0.4, muh = 0.00003), 
          d=list(T=1:T2, p=0),
          trace=TRUE, 
          nls.control(warnOnly=TRUE), 
          algorithm = "port", 
          lower=c(0,0,0,0,0))
          # lower = c(0, 0, 0, 0, 0))
  # d=list() sert pour garder des paramètrse constants 
  # warn.only=TRUE permet de stocker le dernier résultat obtenu même si l'algo échoue
  # l'utilisation de l'algo "port" permet de fixer des bornes à l'espace de recherche des paramètres (ici on veut les garder positifs par exemple)
  # nls2 permet de relancer plusieurs fois l'optimisation en partant de différents points initiaux
  # (note : algos brute-force et grid-search ou random-search font seulement des échantillonages de la fonction objectif)

summary(res2)
par(mfrow=c(1,2))

N = 10001                          # nombre de susceptible
Tmax = 50                  # temps maximal
# par <- list(betah=0.4, betav=0.4, bitting=1, gamma=0.167) # vecteurs de paramètres 
par <- list(muh=coefficients(res2)["muh"], muv=coefficients(res2)["muv"],
            p=0, Cvh=coefficients(res2)["Cvh"], Chv=coefficients(res2)["Chv"], 
            gammah=coefficients(res2)["gammah"])
state <- c(Sh=10000, Ih=1, Rh=0, Sv=50000, Iv=1)              # valeurs initiales du vecteur d'états
comparts2 <- ode(state, 1:Tmax, SIRSI, parms=par)


plot(comparts[,"time"], comparts2[,"Ih"],type="o",xlab="Temps",ylab="Ih(t)", lwd=2, col="blue",  sub="Prédiction dans la population des moustiques")
lines(comparts[,"time"], comparts[,"Ih"],type="l",xlab="Temps",ylab="Ih(t)", col="red", pch = 1)

plot(comparts[,"time"], comparts2[,"Iv"],type="o",xlab="Temps",ylab="Iv(t)", lwd=2, col="blue",  sub="Prédiction dans la population des moustiques")
lines(comparts[,"time"], comparts[,"Iv"],type="l",xlab="Temps",ylab="Iv(t)", col="green", pch = 2)






# # estimation des paramètres p : #BIZZARRE
# res2<-nls2(c(obs, obs2)~model_tofit2(muh, muv, p, Cvh, Chv, gammah, T), 
#           start=list(Chv=0.3, muv = 0.2), 
#           d=list(T=1:Tmax, p=0, muh=1/25000, 
#                     Cvh=coefficients(res)["Cvh"], gammah=coefficients(res)["gammah"]), 
#           trace=TRUE, 
#           nls.control(warnOnly=TRUE), 
#           algorithm = "port", 
#           lower=c(0, 0))
# # d=list() sert pour garder des paramètrse constants 
# # warn.only=TRUE permet de stocker le dernier résultat obtenu même si l'algo échoue
# # l'utilisation de l'algo "port" permet de fixer des bornes à l'espace de recherche des paramètres (ici on veut les garder positifs par exemple)
# # nls2 permet de relancer plusieurs fois l'optimisation en partant de différents points initiaux
# # (note : algos brute-force et grid-search ou random-search font seulement des échantillonages de la fonction objectif)
# 
# summary(res2)
# plot(1:Tmax, obs2, xlab="Temps",ylab="Iv(t)", col="green", pch=4) 
# lines(1:Tmax, fitted(res2),lwd=2, col="blue")
# 
# 
# 
# 
# # estimation des paramètres p :
# res3<-nls2(c(obs, obs2)~model_tofit2(muh, muv, p, Cvh, Chv, gammah, T), 
#            start=list(muh=0.00001), 
#            d=list(T=1:Tmax, p=0,
#                      Cvh=coefficients(res)["Cvh"], gammah=coefficients(res)["gammah"],
#                      Chv=coefficients(res2)["Chv"], muv = coefficients(res2)["muv"]),
#            trace=TRUE, 
#            nls.control(warnOnly=TRUE), 
#            algorithm = "port", 
#            lower=c(0, 0))
# # d=list() sert pour garder des paramètrse constants 
# # warn.only=TRUE permet de stocker le dernier résultat obtenu même si l'algo échoue
# # l'utilisation de l'algo "port" permet de fixer des bornes à l'espace de recherche des paramètres (ici on veut les garder positifs par exemple)
# # nls2 permet de relancer plusieurs fois l'optimisation en partant de différents points initiaux
# # (note : algos brute-force et grid-search ou random-search font seulement des échantillonages de la fonction objectif)
# 
# summary(res3)
# plot(1:Tmax, obs3, xlab="Temps",ylab="Iv(t)", col="green") 
# lines(1:Tmax, fitted(res3),lwd=2, col="blue")





# # facultatif : prédiction et confiance
# confint(res)  # intervalles de confiance sur les estimées des paramètres
# conf<-predict(as.lm(res), interval="confidence")
# lines(d$t, conf[,"lwr"], col="green")
# lines(d$t, conf[,"upr"], col="green")
# conf<-predict(as.lm(res), interval="prediction")
# lines(d$t, conf[,"lwr"], col="dark green")
# lines(d$t, conf[,"upr"], col="dark green")





# ------------------ UTILISATION DES VALEURS ESTIMEES
gammah_estim = coefficients(res)["gammah"]
Cvh_estim = coefficients(res)["Cvh"]
Chv_estim = coefficients(res)["Chv"]
muv_estim = coefficients(res)["muv"]
muh_estim = coefficients(res)["muh"]
# Ieq = N-gamma_estim/beta_estim 
# R0 = beta_estim*N/gamma_estim

  # prédiction à Tmax+1 :
I_pred = model_tofit2(muh=muh_estim, muv=muv_estim, p=0, Cvh=Cvh_estim, Chv=Chv_estim,
                     gammah_estim, 1:(Tmax+1))[Tmax+1]
I_pred




# ------------------ 
# 2°) que se passe-t-il si on n'utilise que les premières données ?
# Par exemple t<= Tmax/2
T2 = round(Tmax/3)
timeline <- 1:T2

obs3 <- comparts[,"Ih"][1:T2] + + rnorm(Tmax, mean=0, sd=2)
obs4 <- comparts[,"Iv"][1:T2] + + rnorm(Tmax, mean=0, sd=2)
res2<-nls2(c(obs3, obs4)~model_tofit2(muh, muv, p, Cvh, Chv, gammah, T), 
          start=list(gammah=0.3, Cvh=0.7, muv = 0.2, Chv=0.35), 
          d=list(T=1:Tmax, p=0, muh = 1/25000),
          trace=TRUE, 
          nls.control(warnOnly=TRUE), 
          algorithm = "port", 
          lower=c(0.2, 0.6, 0.1, 0.3), upper=c(0.4, 0.9, 0.4, 0.5))
# lower = c(0, 0, 0, 0, 0))
# d=list() sert pour garder des paramètrse constants 
# warn.only=TRUE permet de stocker le dernier résultat obtenu même si l'algo échoue
# l'utilisation de l'algo "port" permet de fixer des bornes à l'espace de recherche des paramètres (ici on veut les garder positifs par exemple)
# nls2 permet de relancer plusieurs fois l'optimisation en partant de différents points initiaux
# (note : algos brute-force et grid-search ou random-search font seulement des échantillonages de la fonction objectif)

summary(res2)

par(mfrow=c(1,2))
plot(1:Tmax, obs, xlab="Temps",ylab="Ih(t)", col="red", pch=1) 
lines(1:Tmax, fitted(res2)[1:Tmax],lwd=2, col="blue")

plot(1:Tmax, obs2, xlab="Temps",ylab="Iv(t)", col="red", pch=2) 
lines(1:Tmax, fitted(res2)[(Tmax+1):(2*Tmax)], lwd=2, col="blue")

gammah_estim2 = coefficients(res2)["gammah"]
Cvh_estim2 = coefficients(res2)["Cvh"]
Chv_estim2 = coefficients(res2)["Chv"]
muv_estim2 = coefficients(res2)["muv"]
muh_estim2 = coefficients(res2)["muh"]
# Ieq2 = N-gamma_estim2/beta_estim2 
# R02 = beta_estim2*N/gamma_estim2




# ------------------ 3°) que se passe-t-il si on utilise les données estimées à partir de l'échantillonage ?
Nech=100000
plot(d$t, d$Iech/Nech, ylab = "I/N", col="red", ylim=range(0, 1.5*d$Iech/Nech,na.rm=TRUE ))
points(d$t, d$I/N, col="black")
# intervalle de confiance à 95% :
# conditions : f=N_i/Nech n'est pas proche de 0 ou de 1 ; N_i >= 5 ; Nech-Ni>=5
# (conditions pas toujours vérifiées !!)
q=qnorm(1-0.05/2)
delta = 2*q*sqrt(d$Iech/Nech*(1-d$Iech/Nech)/Nech)
Ilow = d$Iech/Nech - delta/2 
Iup  = d$Iech/Nech + delta/2 
lines(d$t, Ilow, col="red")
lines(d$t, Iup, col="red")
# très larges intervalles de confiance, comme attendu.
# puis refaire calculs avec Iech/Nech *N comme données au lieu de I.
# Ne garder que les valeurs non NA :
obs3 <- d$Iech*N/Nech
timeline3 <- d$t[which(!is.na(obs3))]
obs3 <- obs3[!is.na(obs3)]
res3<-nls(obs3~model_tofit2(beta, gamma, T), start=list(gamma=4, beta=6), d=list(T=timeline3), trace=TRUE, nls.control(warnOnly=TRUE))   
summary(res3)
plot(timeline3, obs3, xlab="Temps",ylab="I(t)", col="red") 
lines(timeline3, fitted(res3),lwd=2, col="blue")
# En enlevant la dernière mesure, peut converger mais erreurs standard énormes !
Ieq = N*(1-coefficients(res3)["gamma"]/coefficients(res3)["beta"])  
R0 = coefficients(res3)["beta"]/coefficients(res3)["gamma"]  

# interprétation biologique : 
# En régime permanent : la proba de guérir correspond à celle de transmettre le virus à un voisin
# mais le problème est qu'il y a la phase d'émergence.
# cette phase fait aussi que la durée moyenne d'infection est biaisée.
# La première personne infectée le reste pendant une durée très longue, 
# qui diminue pour les personnes suivantes.



# ------------------  4°) Analyse de sensibilité globale du modèle
# fonction qui prend en entrée un "d-frame" de paramètres 
# (les X[i,p] où 1<=i<=n où n nombre de paramètres ; 1<=p<=Nrep où Nrep nombre de répétitions)
# et retourne un vecteur Yp de sorties de taille Nrep.
# (une sortie par ligne du d-frame)
par(mfrow=c(1,2))
simulation <- function(paramMatrix){
  V <- c(Sh=N-1, Ih=1, Rh=0, Sv=50000-1, Iv=1)   
  out <-NULL # j'initialise le vecteur renvoyé
  # on parcourt chacune des lignes du d-frame : une ligne = 1 valeur d'un jeux de paramètres.
  for (i in 1:nrow(paramMatrix)){
    # simulation time 
    tps <- seq(1, max(d$t),by=0.5)
    # resolution de l'ODE pour la ième ligne de valeurs de paramètres :
    sol <- ode(V, tps, SIRSI, parms=c("muh"=paramMatrix$muh[i],"muv"=paramMatrix$muh[i],
              "Cvh"=paramMatrix$Cvh[i], "Chv"=paramMatrix$Chv[i],"gammah"=paramMatrix$gammah[i],
              "p"=paramMatrix$p[i]), method="lsoda")
    # Remarque : on peut aussi utiliser la méthode rk4 pour gagner en précision, si besoin. Voir l'aide de 'ode'.
    # ici j'ai choisi comme variable d'intérêt le nombre d'infectés au temps final
    # mais on pourrait imaginer d'autres variables comme par exemple la somme des carrés
    # utilisées pour model_tofit pour l'estimation :
    out<-c(out, sol[length(tps),"Ih"])
  }
  return(out)
}

Nrep=10000 # passer à plutot 10 000 une fois que c'est testé, pour avoir convergence.
# vérifier la convergence en lançant plusieurs fois la méthode Sobol et regardant si résultats similaires. 

# Matrice d'échantillonage :
X1<-data.frame(muv=runif(Nrep, 0, 6), muh=runif(Nrep, 0, 6), gammah=runif(Nrep, 0, 6),
               Chv=runif(Nrep, 0, 6), Cvh=runif(Nrep, 0, 6), p=runif(Nrep, 0, 6))
# Matrice de ré-échantillonage :
X2<-data.frame(muv=runif(Nrep, 0, 6), muh=runif(Nrep, 0, 6), gammah=runif(Nrep, 0, 6),
               Chv=runif(Nrep, 0, 6), Cvh=runif(Nrep, 0, 6), p=runif(Nrep, 0, 6))
# remarque: on peut pour certains paramètres mettre une valeur constante, 
# ce qui revient à considérer qu'ils sont connus avec précision et 
# donc leur incertitude n'est pas prise en compte dans le calcul des indices.

# Faire un petit test de la fonction simulation, avant de passer à la suite.
simulation(X1)

# indices d'ordre 1 et total
res2<-sobol2007(simulation, X1, X2, nboot=30)
# res2<-sobol2002(simulation, X1, X2, nboot=10)
res18 <- res2
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
# La fonction définissant le modèle pour Sobol prend en argument une matrice de taille Nrep x p dont chaque ligne représente une réalisation 
# d'un jeu de p paramètres, tirée selon les lois définies pour caractériser l'incertitude de chacun d'entre eux. A l'intérieur de cette fonction,
# on va lancer autant de simulations qu'on donne de jeux de paramètres différents, c'est-à-dire Nrep, le nombre de lignes de la matrice (nrow()).
# Pour chaque simulation, il faut définir une sortie scalaire (ie ce doit être un réel : si on en veut plusieurs, il faut lancer autant de fois l'analyse). 
# Dans mon exemple, j'ai choisi de prendre la valeur du nombre d'infectés au dernier temps de ma simulation (Tmax). A vous de choisir une sortie d'intérêt
# pour votre modèle. La fonction doit renvoyer un vecteur dont chaque composante correspond à la sortie d'une simulation : ce vecteur est donc de taille Nrep.
# Ensuite la fonction Sobol2002/7 prend en argument le nom de cette fonction, et deux matrices de jeux de paramètres X1 et X2. Elle comprend aussi un argument de
# nombre de répétition de "bootstrap" qui sert à estimer l'intervalle de confiance sur la valeur de l'indice.

# Pour méthode Morris : faire un print la 1ère fois permet de voir la
# matrice d'exploration de l'espace des paramètres (l'enlever ensuite).
simulationM <- function(paramMatrix){
  V <- c(Sh=10000-1, Ih=1, Rh=0, Sv=50000-1, Iv=1)
  out <-NULL
  
  for (i in 1:nrow(paramMatrix)){
    # parameters
    # ici param doit être au format liste (indifférent dans méthode précédente)
    par <- list("muh"=paramMatrix[i,1],"muv"=paramMatrix[i,2],"p"=paramMatrix[i,3],
                "Cvh"=paramMatrix[i,4], "Chv"=paramMatrix[i,5],"gammah"=paramMatrix[i,6]) 
    #print(unlist(par))
    tps <- seq(1, max(d$t),by=0.1)
    sol <- ode(V, tps, SIRSI, par, method="lsoda") 
    out<-c(out, sol[length(tps),"Ih"])
  }
  return(out)
}
resM <- morris(model = simulationM, factors = 6, r = 50, binf=0, bsup=1,
               design = list(type = "oat", levels = 10, grid.jump = 1))
plot(resM, sub="Analyse de sensibilité globale pour le modèle SIRSI selon la méthode de Morris")
resM



