# This version of the algorithm is the better way of doing it where the distances
# of the two changed clusters only are calculated each time we make an iteration.

#This function takes a vector v and scales it so the max value gets a 1 and
#the lowest a 0. It gives back a scaled vector


scale.01 <- function(v){
  sc.01 <- (v-min(v))/(max(v)-min(v))
  sc.01
}


# This function calculates the energy of a single cluster only (coi- cluster of interest)
calc.Ek <- function(m,clus,coi){

  clus.d <- m[which(clus==coi),]
  Ek <- sum(dist(clus.d))
  Ek
}

# We now use the function above to make one that calculates the E for each cluster
calc.E.all <- function(m,clus){
  
  Eks <- vector("numeric",max(clus)) # vactor to catch the energies
  
  for(i in 1:max(clus)){
    
    Eks[i] <- calc.Ek(m,clus,i) # calc the energy of cluster i
    
  }
  Eks # return the energies
  
}

# a simple function to cal;culate the mean energy according to the formulas.
# mean() would be the same
E.tot <- function(Eks){
  
  sum(Eks)/length(Eks)
  
}

# This function takes Eold and Enew and calculated the exponent to be compared
# against a random number.
calc.exp <- function(E.new,E.old,T){
  exp(-((E.new-E.old)/T))
}


#This is the main algorithm that performs the annealing. It takes the data,how
#many K we are looking for, The number of iterations to perform, starting
#temperature and the cooling factor.

sa.ok <- function(data,K,Iter,Temp,cool){

  clusters <- sample(1:K,nrow(data),replace = T) #initialise random clusters
  
  #clusters.o <- clusters
  
  Es.old <- calc.E.all(data,clusters)
  #E.old <- E.tot(Es.old)
 
  for(i in 1:Iter){   # start iterating 
    
    clusters.new <- clusters #copy the clusters
    #Es.new <- Es.old
    
    row.id <- sample(1:nrow(data),1) #pick a gene at random
    
    from.c <- clusters.new[row.id] # get the cluster it's moving from
    #to.c <- sample((1:K)[!(1:K) %in% from.c],1)
    to.c <- sample((1:K)[-from.c],1) # randomly choose a new cluster
    
    clusters.new[row.id] <-  to.c # replace the old cluster with the new
    
    Es.new <- Es.old #make a copy of the energies vector
    # calc the energies of the two changed clusters
    Es.new[from.c] <- calc.Ek(data,clusters.new,from.c) 
    Es.new[to.c] <- calc.Ek(data,clusters.new,to.c)
    
    E.new <- E.tot(Es.new) # calculate the new average E
    E.old <- E.tot(Es.old) # calculate the old average E
    
    if(E.new < E.old){  # if new < old accept the move copy the new clusters into the previous
      clusters <-  clusters.new #copy the new clusters into the previous
      Es.old <- Es.new # make Enew to Eold
    }else{
      
      if(calc.exp(E.new,E.old,Temp) > runif(1)){ #evaluate the exprssion against the random number from runif(1)
        clusters <- clusters.new #copy the new clusters into the previous
        Es.old <- Es.new  # make Enew to Eold
      }
    }
    
    {cat("\r",E.old)} #print out the energy to the screen
    
    Temp <- Temp*cool # cool the system
  }
  clusters # return the clusters
}

ycc <- read.delim("../Spellman_Yeast_Cell_Cycle.tsv",row.names=1,header=T,sep="\t")
ycc <- as.matrix(ycc)
ycc.01 <- t(apply(ycc,1,scale.01))




system.time(proc.time(clus <- sa.ok(ycc.01,10,25000,20,0.995)))

# plot the clusters using base R
par(mfrow=c(3,4))

for(i in 1:max(clus)){
  
  ycc.c <- ycc.01[which(clus==i),]
  plot(ycc.c[1,],ty="l",ylim=range(ycc.c))
  apply(ycc.c,1,lines)
  
}

