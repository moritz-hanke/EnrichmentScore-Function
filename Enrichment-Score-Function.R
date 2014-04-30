# Random p-Werte eines Chi^2-Verteilung erstellen
set.seed(123)
p.zahlen.rand <- matrix(pchisq(rchisq(20000, 2,2), df=2, ncp=2), ncol=20, byrow=T)
set.seed(321)
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



max_ES <- function(set, variables, daten){   ### WICHTIG: Spalten müssen SNPs sein, Zeilen die
                                               # beobachteten bzw. permutierten Teststatistiken
                                             ### Observed data have to be a 1row, m-column matrix
  
  if((max(daten) > 1) | (min(daten) < 0)){
    stop("Are you realy using p-values? Values don't lie between 0 and 1")
  }
  
  if((is.matrix(daten) & is.matrix(variables)) ==FALSE){
    stop("Observed data have to be a Matrix with one row and a column for every SNP; 
         Permutated data have to be a Matrix with rows=number_of_permutations and
         a column for every SNP")
  }
  
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
  ES_norm <- ES_max/length(set)    # maxES wird mit Groesse des PW normalisiert
  ES_norm
}

#var <- zahlen.rand[1,]
#var[names(var) %in% set1]   # werte eines sets
#seq_along(var)[names(sort(var, decreasing=T)) %in% set1] # rausfinden, welche positionen ein set in der sortierung belegt
#var[names(sort(var, decreasing=T)) %in% set1] # irgendwas verkehrt
#sort(var, decreasing=T)[names(sort(var, decreasing=T)) %in% set1]




# für set3
obs.manipuliert1 <- apply(p.zahlen.obs.manipuliert, MARGIN=1, max_ES, set=set1, daten=p.zahlen.obs.manipuliert) #wichtig:
                                                                # X von apply muss daten entsprechen und
                                                                # MARGIN =1 sein
obs.manipuliert2 <- apply(p.zahlen.obs.manipuliert, MARGIN=1, max_ES, set=set2, daten=p.zahlen.obs.manipuliert)
obs.manipuliert3 <- apply(p.zahlen.obs.manipuliert, MARGIN=1, max_ES, set=set3, daten=p.zahlen.obs.manipuliert)


obs1 <- apply(p.zahlen.obs, MARGIN=1, FUN=max_ES, set=set1, daten=p.zahlen.obs)
obs2 <- apply(p.zahlen.obs, MARGIN=1, FUN=max_ES, set=set2, daten=p.zahlen.obs)
obs3 <- apply(p.zahlen.obs, MARGIN=1, FUN=max_ES, set=set3, daten=p.zahlen.obs)

rand1 <- apply(p.zahlen.rand, MARGIN=1, max_ES, set=set1, daten=p.zahlen.rand)
rand2 <- apply(p.zahlen.rand, MARGIN=1, max_ES, set=set2, daten=p.zahlen.rand)
rand3 <- apply(p.zahlen.rand, MARGIN=1, max_ES, set=set3, daten=p.zahlen.rand)




sum((obs1 < rand1))/length(rand1)    ### alle rand sind gleich lang -> Anzahl Permutationen
sum((obs2 < rand2)/length(rand2))
sum((obs3 < rand3)/length(rand3))

sum((obs.manipuliert1 < rand1)/length(rand1))
sum((obs.manipuliert2 < rand2)/length(rand2))
sum((obs.manipuliert3 < rand3)/length(rand3))


### Sets müssen als Liste vorliegen; Namen der Elemente der Liste
  # sollten PW/Sets entsprechen!
Sets <- list(set1,set2,set3)
names(Sets) <- unlist(strsplit("set1,set2,set3", ","))


ES <- function(data.obs, data.perm, sets){
  ### Fehler abfangen
  if(is.null(names(Sets))){
    stop("List which contains the Pathways/Sets should have a name for each element")
  }
  if(is.data.frame(data.obs) | is.matrix(data.obs) != TRUE){      
    stop("data.obs must be dataframe or matrix")
  }
  if(is.data.frame(data.perm) | is.matrix(data.perm) != TRUE){
    stop("data.perm must be dataframe or matrix")
  }
  if(is.list(sets)!=TRUE){
    stop("sets has to be a list")
  }
  if(length(ncol(data.obs))!=length(ncol(data.perm))){
    stop("Length of data.ob and data.perm differ")
  }
  if(!all(sort(colnames(data.obs)) == sort(colnames(data.perm)))){
    stop("data.obs and data.perm have different colnames")
  }
  ### Fehler abfangen Ende
  
  results <- vector(mode="list", length=length(sets))
  
  maxES.obs <- rep(NA, length(sets))    ### Erstellen von Objekten, die die maxEX der
  names(maxES.obs) <- sets                # beobachteten und permutierten Daten enthalten
                                          # sollen
  maxES.perm <- vector(mode="list", length=length(sets))  
  names(maxES.perm) <- sets
  
  maxES.obs <- sapply(sets, FUN=function(x){   ### für alle Sets/PW soll Max_ES bei den
                                                 # beobachteten Daten angewedet werden
    apply(data.obs, MARGIN=1, FUN=max_ES, set=x, daten=data.obs)
  }
  )
  
  maxES.perm <- sapply(sets, FUN=function(x){   ### für alle Sets/PW soll Max_ES bei den
                                                  # permutierten Daten angewedet werden
    apply(data.perm, MARGIN=1, FUN=max_ES, set=x, daten=data.perm)
  }
  )
  
  
  nperm <- length(data.perm[,1]) ### wieviele Permutationen gab es? fuer normierten P-Wert
  
  nom_P <- sapply(colnames(maxES.perm), FUN=function(X){ ### Berechnung des Normierten
                                                           # P-Werts
    (sum(maxES.obs[X] < maxES.perm[,X]))/nperm
  })
  back <- list(nom_P, maxES.obs) ### Rueckgabe der normierten P-Werte
  names(back) <- c("nominal p-value", "normalized max_ES")
  back
}

system.time(
ergebnisse <- ES(data.obs=p.zahlen.obs, data.perm=p.zahlen.rand, Sets)
)

