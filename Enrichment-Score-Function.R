# Random zahlen erstellen

zahlen.rand <- matrix(abs(round(rchisq(2000, 2,2), 2)), ncol=20, byrow=T)
zahlen.obs <- matrix(abs(round(rchisq(20, 2,2), 2)), ncol=20, byrow=T)
colnames(zahlen.rand) <- letters[1:20]
colnames(zahlen.obs) <- letters[1:20]

zahlen.obs.manipuliert <- zahlen.obs
zahlen.obs.manipuliert[,14:20] <- c(0.03,0.01,0.03,0.07,0.01,0.02, 0.01)
colnames(zahlen.obs.manipuliert) <- letters[1:20]

### Sets sind in diesem Fallnicht ein Gen sondern ein Pathway, der SNPs enthält
set1 <- letters[c(1:4, 6)]
set2 <- letters[c(5,7:13)]
set3 <- letters[14:20]


ES <- function(var, set, daten){   # wichtig: spalten müssen SNPs sein, Zeilen die beobachteten
                                                    # bzw. permutierten teststatistiken
  temp <- var
  names(temp) <- names(daten)
  temp <- sort(var, decreasing=T)
  N <- ncol(daten)
  N_h <- length(set)
  N_r <- sum(var[names(var) %in% set]) # summe der hits eines Sets
  
  Hits <- rep(0, length(temp)) 
  Hits[seq_along(var)[names(sort(var, decreasing=T)) %in% set]] <- ### Vektor,der nur 0 und an den passenden
    sort(var[names(var) %in% set], decreasing=T)                     # sortierten Positionen die Werte von 
                                                                     # einem Set enthält 
  Hits <- Hits/N_r                                                  
  
  Miss <- rep( 1/(N-N_h), N)
  Miss[seq_along(var)[names(sort(var, decreasing=T)) %in% set]] <- 0
  
  Phit <- cumsum(Hits)
  Pmiss <- cumsum(Miss)
  
  ES <- Phit - Pmiss
  ES_max <-max(ES)
  
}

#var <- zahlen.rand[1,]
#var[names(var) %in% set1]   # werte eines sets
#seq_along(var)[names(sort(var, decreasing=T)) %in% set1] # rausfinden, welche positionen ein set in der sortierung belegt
#var[names(sort(var, decreasing=T)) %in% set1] # irgendwas verkehrt
#sort(var, decreasing=T)[names(sort(var, decreasing=T)) %in% set1]


# für set3
obs.manipuliert <- apply(zahlen.obs.manipuliert, MARGIN=1, ES, set=set3, daten=zahlen.obs.manipuliert) #wichtig:
                                                                # X von apply muss daten entsprechen und
                                                                # MARGIN =1 sein

obs <- apply(zahlen.obs, MARGIN=1, ES, set=set3, daten=zahlen.obs)

rand <- apply(zahlen.rand, MARGIN=1, ES, set=set3, daten=zahlen.rand)



