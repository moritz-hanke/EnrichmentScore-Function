# Random zahlen erstellen
zahlen.rand <- matrix(rchisq(2000, 2,2), ncol=20, byrow=T)
zahlen.obs <- matrix(rchisq(20, 2,2), ncol=20, byrow=T)
colnames(zahlen.rand) <- letters[1:20]
colnames(zahlen.obs) <- letters[1:20]

zahlen.obs.manipuliert <- zahlen.obs
zahlen.obs.manipuliert[,14:20] <- c(0.03,0.01,0.03,0.07,0.01,0.02, 0.01)
colnames(zahlen.obs.manipuliert) <- letters[1:20]

# Random p-Werte erstellen
p.zahlen.rand <- matrix(pchisq(rchisq(2000, 2,2), df=2, ncp=2), ncol=20, byrow=T)
p.zahlen.obs <- matrix(pchisq(rchisq(20, 2,2), df=2, ncp=2), ncol=20, byrow=T)
colnames(p.zahlen.rand) <- letters[1:20]
colnames(p.zahlen.obs) <- letters[1:20]

p.zahlen.obs.manipuliert <- p.zahlen.obs
p.zahlen.obs.manipuliert[,14:20] <- c(0.9987,0.98778,0.9934,0.96434,0.9874,0.82345, 0.99)
p.zahlen.obs.manipuliert[,14:20] <- 1-c(0.9987,0.98778,0.9934,0.96434,0.9874,0.9982345, 0.99)
colnames(p.zahlen.obs.manipuliert) <- letters[1:20]


### Sets sind in diesem Fallnicht ein Gen sondern ein Pathway, der SNPs enthält
set1 <- letters[c(1:4, 6)]
set2 <- letters[c(5,7:13)]
set3 <- letters[14:20]


####################
### QUANTILFUNKTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
####################
### QUANTILFUNKTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
####################
### QUANTILFUNKTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
####################
### QUANTILFUNKTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!



ES <- function(variables, set, daten){   # wichtig: spalten müssen SNPs sein, Zeilen die beobachteten
                                                    # bzw. permutierten teststatistiken
  temp <- qnorm(1-variables)            ### z-Scores der p-werte mittels Quantilfunktion der SNV
  namen <- names(variables)[order(temp, decreasing=T)] ### namen der nachfolgend sortierten
  # z-Scores
  temp <- sort(temp, decreasing=T)  ### Sortierung der z-Scores
  names(temp) <- namen                # Namesgebung der z-Scores
  
  N <- ncol(daten)           ### wie viele SNPs/Spalten gibt es?
  N_h <- length(set)         ### Wieviele SNPs enthält das Set/der Pathway?
  N_r <- sum(abs(temp[names(temp) %in% set])) ### Summe der z-Scores aus einem Sets
  
  Hits <- rep(0, length(temp)) 
  Hits[seq_along(temp)[names(temp) %in% set]] <- ### Vektor,der nur 0 und an den passenden
    abs(temp[names(temp) %in% set])                # sortierten Positionen die Werte von 
  # einem Set enthält 
  Hits <- Hits/N_r        ### mit N_r normierten sortierten Werte
  
  Miss <- rep( 1/(N-N_h), N)  ### Gegenstrück zu Hits; alle Einträge die != 0 sind gleich
  Miss[seq_along(temp)[names(temp) %in% set]] <- 0
  
  Phit <- cumsum(Hits)
  Pmiss <- cumsum(Miss)
  
  ES_all <- Phit - Pmiss
  ES_max <-max(ES_all)
  ES_max
  
}

#var <- zahlen.rand[1,]
#var[names(var) %in% set1]   # werte eines sets
#seq_along(var)[names(sort(var, decreasing=T)) %in% set1] # rausfinden, welche positionen ein set in der sortierung belegt
#var[names(sort(var, decreasing=T)) %in% set1] # irgendwas verkehrt
#sort(var, decreasing=T)[names(sort(var, decreasing=T)) %in% set1]




# für set3
obs.manipuliert <- apply(p.zahlen.obs.manipuliert, MARGIN=1, ES, set=set1, daten=p.zahlen.obs.manipuliert) #wichtig:
                                                                # X von apply muss daten entsprechen und
                                                                # MARGIN =1 sein

obs1 <- apply(p.zahlen.obs, MARGIN=1, FUN=ES, set=set1, daten=p.zahlen.obs)
obs2 <- apply(p.zahlen.obs, MARGIN=1, FUN=ES, set=set2, daten=p.zahlen.obs)
obs3 <- apply(p.zahlen.obs, MARGIN=1, FUN=ES, set=set3, daten=p.zahlen.obs)

rand1 <- apply(p.zahlen.rand, MARGIN=1, ES, set=set1, daten=p.zahlen.rand)
rand2 <- apply(p.zahlen.rand, MARGIN=1, ES, set=set2, daten=p.zahlen.rand)
rand3 <- apply(p.zahlen.rand, MARGIN=1, ES, set=set3, daten=p.zahlen.rand)

plot()


sum((obs < rand1)/length(rand1))    ### alle rand sind gleich lang -> Anzahl Permutationen
sum((obs < rand2)/length(rand2))
sum((obs < rand3)/length(rand3))

NES.obs1 <- obs/mean(rand1)
NES.obs2 <- obs/mean(rand2)
NES.obs3 <- obs/mean(rand3)
NES.rand1 <- rand1/mean(rand1)
NES.rand2 <- rand2/mean(rand2)
NES.rand3 <- rand3/mean(rand3)

hist(rand)

sum(NES.obs1 < NES.rand1)/length(NES.rand1)       ### kleiner < weil man den Anteil wissen will,
                                                   # der durch den Zufall einen größeren ES hat
sum(NES.obs2 < NES.rand2)/length(NES.rand2)
sum(NES.obs3 < NES.rand3)/length(NES.rand3)

