# random zahlen erstellen

zahlen.rand <- matrix(abs(round(rnorm(2000, 4, 3), 2)), ncol=20, byrow=T)
zahlen.obs <- matrix(abs(round(rnorm(20, 4, 3), 2)), ncol=20, byrow=T)
colnames(zahlen.rand) <- letters[1:20]
colnames(zahlen.obs) <- letters[1:20]

zahlen.obs.manipuliert <- zahlen.obs
zahlen.obs.manipuliert[1,14:19] <- c(0.03,0.01,0.03,0.07,0.01,0.02)
colnames(zahlen.obs.manipuliert) <- letters[1:20]
set1 <- letters[c(1:4, 6)]
set2 <- letters[c(5,7:13)]
set3 <- letters[14:20]

# Frage: kann es eigentlich, wenn man Chi^2 Werte nimmt, eine Ansammlung am Ende geben?
#Ich meine ja, wenn sie alle sehr klein sind. Das Problem dabei ist, dass das dann ja
# eigentlich sehr hohen P-Werten entspricht und folglich würde man nicht von einer
# Assoziation ausgehen
# können sich aber logisch die werte nur an einer seite (oben, weil große werte) sammeln,
# stören sie sich so ggf in ihrem signal: bspw. wäre als genexpression eine starke positive
# und eine starke negative korrelation zwischen zwei pathways und dem phaenotypen zu sehen,
# aber durch den chi^2 sind es nur noch starke "positive" "korrelationen"

ES <- function(var, set, daten){   # wichtig: spalten müssen SNPs sein, Zeilen die beobachteten
                                                    # bzw. permutierten teststatistiken
  temp <- var
  names(temp) <- names(daten)
  temp <- sort(var, decreasing=T)
  N <- ncol(daten)
  N_h <- length(set)
  N_r <- sum(var[names(var) %in% set]) # summe der hits eines Sets
  
  Hits <- rep(0, length(temp)) 
  Hits[seq_along(var)[names(sort(var, decreasing=T)) %in% set]] <- sort(var[names(var) %in% set], decreasing=T)# Vektor,der nur 0 und an den passenden sortierten Positionen die Werte von einem Set enthält
  Hits <- Hits/N_r
  
  Miss <- rep( 1/(N-N_h), N)
  Miss[seq_along(var)[names(sort(var, decreasing=T)) %in% set]] <- 0
  
  Phit <- cumsum(Hits)
  Pmiss <- cumsum(Miss)
  
  ES <- Phit - Pmiss
  ES
  #ES_max <- max(ES)
  #ES_max
  
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



names(zahlen.rand[4,5])
