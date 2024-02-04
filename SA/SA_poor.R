
# This version of the algorithm is the bad way of doing it where the distances
#of each cluster are calculated each time we make an iteration.

#This function takes a vector v and scales it so the max value gets a 1 and
#the lowest a 0. It gives back a scaled vector

scale.01 <- function(v){
    sc.01 <- (v-min(v))/(max(v)-min(v))
    sc.01
}

#This function does everything. It takes a vector of cluster assignments and
#loops through each one in i. The data is then reduced to just thos genes 
#belonging to clusters i, and then the poirwise distances are calculated and
#summed. The total E is then ccalculated.

calc.E.tot <- function(m,clus){

    clus.E <- NULL # an empty vector to catch the energies

    for(i in 1:max(clus)){ #loop through each clusters
        #print(i)
        clus.d <- m[which(clus==i),] #reduce to those belonging to cluster i
        #clus.V <- c(clus.V,sqrt(sum(dist(clus.d)^2)))
        clus.E <- c(clus.E,sum(dist(clus.d))) #calc distances and sum. put E in clus.V
    }
    
    sum(clus.E)/max(clus) # calculate the average energy and return it.
}


# This function takes Eold and Enew and calculated the exponent to be compared
#against a random number.

calc.exp <- function(E.new,E.old,T){
    exp(-((E.new-E.old)/T))
}

#This is the main algorithm that performs the annealing. It takes the data,how
#many K we are looking for, The number of iterations to perform, starting
#temperature and the cooling factor.


sa.nsg <- function(data,K,Iter,Temp,cool){
  
  
  clusters <- sample(1:K,nrow(data),replace = T) #initialise random clusters
  
  E.old <- calc.E.tot(data,clusters) # calculate the energy
  
  for(i in 1:Iter){   # start iterating
  
      clusters.new <- clusters # make a copy pof the current clusters
      #E.old <- calc.V.tot(data,clusters) # calculate the energy
      
  
      row.id <- sample(1:nrow(data),1) #pick a gene at random
      
      clusters.new[row.id] <- sample(1:K,1) # choose a random cluster and place it in the new cluster config
  
      E.new <- calc.E.tot(data,clusters.new) #calculate a new energy
  
      if(E.new==E.old){ # if Enew is the same, cool and sip to the next iteration
          Temp <- Temp*cool #cool the system
          next # go to the next iteration
      }
  
      if(E.new < E.old){ # if new < old accept the move copy the new clusters into the previous
          clusters <-  clusters.new #copy the new clusters into the previous
          E.old <- E.new # make Enew to Eold
          
      }else{
         if(calc.exp(E.new,E.old,Temp) > runif(1)){ #evaluate the exprssion against the random number from runif(1)
          clusters <- clusters.new #copy the new clusters into the previous
          E.old <- E.new  # make Enew to Eold
         }
      }
      {cat("\r",E.old)} #print out the energy to the screen
  
      Temp <- Temp*cool # cool the system
  }
  print(E.old) # print the final Energy
  clusters # return the clusters
  
}

#Lets use these functions

# read the data
ycc <- read.delim("../Spellman_Yeast_Cell_Cycle.tsv",row.names=1,header=T,sep="\t")

#make it a matrix
ycc <- as.matrix(ycc)

#scale the genes using the apply function with out function
ycc.01 <- t(apply(ycc,1,scale.01))

# run the algorithm and measure the time it takes
system.time(clus <- sa.nsg(ycc.01,10,25000,20,0.995))

# compare to a kmeans clustering again
calc.E.tot(ycc.01,kmeans(ycc.01,10)$cluster)

# DID YOU BEAT IT?


# plot the clusters using base R
par(mfrow=c(3,4))

for(i in 1:max(clus)){
  
  ycc.c <- ycc.01[which(clus==i),]
  plot(ycc.c[1,],ty="l",ylim=range(ycc.c))
  apply(ycc.c,1,lines)
  
}

