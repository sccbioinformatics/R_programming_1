# This version of the algorithm is an elegant solution that Elena impemented where
# all the distances between all genes are calculated up front, and then we simply
# subset the distances into clusters and takle the sum. This means there is no 
# need to caluculate distances during the iterations.

# The dist() function returns a diagonal matrix which we then make square using
# as.matrix. The means we can then subset any distance matrix dm using dm[ind,ind]

#This function takes a vector v and scales it so the max value gets a 1 and
#the lowest a 0. It gives back a scaled vector

scale.01 <- function(v){
  sc.01 <- (v-min(v))/(max(v)-min(v))
  sc.01
}

# This function calulates the energy of a single cluster and subsets a square
# distance  matrix dm and returns the sum of it
calc.Ek.lu <- function(dm,clusters,coi){
  
  rind <- which(clusters==coi)
  
  Ek <- sum(dm[rind,rind])/2
  Ek
}

# Takes the function above and does the same fao all clusters storing the enegoies
# in Es

calc.E.lu.all <- function(dm,clus){
  
  Es <- vector("numeric",length=max(clus))
  
  for(i in 1:max(clus)){
    
    Es[i] <- calc.Ek.lu(dm,clus,i)
    
  }
  Es
}

# a simple function to cal;culate the mean energy according to the formulas.
# mean(Eks) would be the same
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
sa.eb <- function(data,K,Iter,Temp,cool){
  
  #caculates the pair-wise distances of all genes and makes a matrix
  data.dst <- as.matrix(dist(data)) 
  
  
  clusters <- sample(1:K,nrow(data),replace = T) #initialise random clusters
  
  Es.old <- calc.E.lu.all(data.dst,clusters) # calc e for each cluster
  
  for(i in 1:Iter){   
    
    clusters.new <- clusters #copy the clusters
    
    row.id <- sample(1:length(clusters),1) #pick a gene at random
    
    from.c <- clusters.new[row.id] # get the cluster it's moving from
    to.c <- sample((1:K)[-from.c],1) # randomly choose a new cluster
    
    clusters.new[row.id] <-  to.c # replace the old cluster with the new
    
    Es.new <- Es.old #make a copy of the energies vector
    
    # calc the energies of the two changed clusters
    Es.new[from.c] <- calc.Ek.lu(data.dst,clusters.new,from.c)
    Es.new[to.c] <- calc.Ek.lu(data.dst,clusters.new,to.c)
    
    E.new <- E.tot(Es.new)
    E.old <- E.tot(Es.old)
    
    if(E.new < E.old){ # if new < old accept the move copy the new clusters into the previous
      clusters <-  clusters.new #copy the new clusters into the previous
      Es.old <- Es.new # make Enew to Eold 
    }else{
      
      if(calc.exp(E.new,E.old,Temp) > runif(1)){ #evaluate the exprssion against the random number from runif(1)
        clusters <- clusters.new #copy the new clusters into the previous
        Es.old <- Es.new  # make Enew to Eold
      }else{
        next
      }
    }
    
    {cat("\r",E.old)} #print out the energy to the screen
    
    Temp <- Temp*cool
  }
  clusters
}

ycc <- read.delim("../Spellman_Yeast_Cell_Cycle.tsv",row.names=1,header=T,sep="\t")
ycc <- as.matrix(ycc)


ycc.01 <- t(apply(ycc,1,scale.01))

system.time(clus <- sa.eb(ycc.01,8,25000,20,0.995))

par(mfrow=c(2,4))

for(i in 1:max(clus)){
  
  ycc.c <- ycc.01[which(clus==i),]
  plot(ycc.c[1,],ty="l",ylim=range(ycc.c))
  apply(ycc.c,1,lines)
  
}

#########
e.dist <- function(v1, v2){
  d <- sqrt(sum((v1-v2)^2))
  d
}

my.kmeans <- function(data,nclus,iter){ # setup a functon that will take some data, number of clusters, and max number of iterations.
  
  centroids <- data[sample(1:nrow(data), nclus, replace=F),] 
  
  dists <- matrix(0, nrow = nrow(data), ncol = nclus) 
  
  
  clusters <- sample(1:nrow(centroids),nrow(data),replace=T) #make a random assignment of clusters to start off
  
  for(iteration in 1:iter){ # start the iterations
    
    for(gene in 1:nrow(data)){
      
      for(k in 1:nclus){
        dists[gene, k] <- e.dist(data[gene,], centroids[k,]) # for each gene calculate the distance to each centroid (k).
      }
    }
    
    clusters.new <- apply(dists, 1, which.min) # get the mnew clusters assignments
    if(all.equal(clusters,clusters.new)==T){ #see if it equals the previous round
      break # if so, stop
    }else(clusters <- clusters.new) # otherwise set the clusters to the new.clusters
    
    for(k in 1:nrow(centroids)){
      centroids[k,] <- apply(data[which(clusters == k),], 2, mean) # define new centroids.
    }
    print(iteration)
  }
  clusters.new # return the clusters.
}

clus.8 <- my.kmeans(ycc.01,8,100) # run it using 5 clusters


