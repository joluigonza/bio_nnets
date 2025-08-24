
library(ggplot2)


#common gems for individuals wrt to a treatment
  
common_gens <- function(x, i_stds, ngens, treat, par, tol) {
  
  i=1
  
  thisx <- x[i_stds[1:ngens],3]
  
  thisy <- thisx-x[i_stds[1:ngens],3-treat]
  #print(length(thisy[(-1)^par*thisy>0]))
  
  thisins <- which((-1)^par*thisy > tol*thisx)
  
  i= i+1
  
  for (j in seq(6,15, by=3)){
    
    thisx <- x[i_stds[1:ngens],j]
    thisy <- thisx-x[i_stds[1:ngens],j-treat]
    
    #print(length(thisy[(-1)^par*thisy>0]))
    
    newins <- which((-1)^par*thisy > tol*thisx)
    
    thisins <- thisins[which(thisins %in% newins)]
    
    i= i+1
  }
  return (i_stds[thisins])
}

##################################################################

null_gens <- function(x, i_stds, ngens, treat, tol) {  
  
  i=1
  
  thisx <- x[i_stds[1:ngens],3]
  
  thisy <- abs(thisx-x[i_stds[1:ngens],3-treat])
  
  #print(length(thisy[(-1)^par*thisy>0]))
  
  thisins <- which( thisy < tol*thisx)
  
  i= i+1
  
  for (j in seq(6,15, by=3)){
    
    thisx <- x[i_stds[1:ngens],j]
    
    thisy <- abs(thisx-x[i_stds[1:ngens],j-treat])
    
    #print(length(thisy[(-1)^par*thisy>0]))
    
    newins <- which(thisy < tol*thisx)
    
    thisins <- thisins[which(thisins %in% newins)]
    
    i= i+1
  }
  return (i_stds[thisins])
}


######################## intersection from a list of genes, very useful 

intersect_gens <- function(list, vec) {
  
  thislist <- list[[vec[1]]]
  for (i in vec[2:length(vec)]) {
    
    thislist <- thislist[which(thislist %in% list[[i]])]
    
  }
  
  return (thislist)
  
}
  
#########################################

build_nn <- function(trainin) {
  
  trainy=matrix(0, ntrain, 1)
  
  i=1
  
  for (j in (1:ntrain)*3 - treat) {
    trainy[i]=x[i_stds[trainin],j]
    i = i+1
  }
  
  testy=matrix(0,ntest,1)
  i=1
  for (j in ((ntrain+1):(ntrain + ntest))*3 - treat) {
    testy[i]=x[i_stds[trainin],j]
    i = i+1
  }
  
  
  
  nn <- nnet(traindx, trainy, size=ng, maxit=10000, MaxNWts = 100000,
       linout=TRUE, act.fct=relu)
  
  
  thiserr <-  max(abs(testy-predict(nn, testdx)))/max(abs(testy))

  while (thiserr > tol) {
    nn <- nnet(traindx, trainy, size=ng, maxit=10000, MaxNWts = 100000,
               linout=TRUE)

    thiserr <-  max(abs(testy-predict(nn, testdx)))/max(abs(testy))

  }

  
  return(nn)
  
  #return(nn$wts)
  
  #return(nn$fitted.values)
  
  
}

###############################################



predict_nn <- function(thisin) {
  
  this <- predict(nns[[thisin]], input)
  thismin <- min(this)
  thismax <- max(this)
  
  if (thismin > min(treat[1:nrow(treat), thisin]) && thismax < max(treat[1:nrow(treat), thisin])) {
    ans <- 1
    
  } else {
    
    ans <- 0
  }
  
  return (ans)
  
}
###################################################



predict_nn1 <- function(thisin) {
  
  this <- predict(nns[[thisin]], input)
  
  return (this)
  
}

####################################################

gene_ratios <- function(x,inds,treat){
  
  thisng = length(inds)
  dy=matrix(0,thisng,1);
  dx=matrix(0,thisng,1);
  
iter=1
for (j in inds) {
  dy[iter]=sort(x[j,(1:5)*3 - treat],decreasing = TRUE)[2]-sort(x[j,(1:5)*3 - treat])[2]
  iter = iter +1
  
}

iter=1
for (j in inds) {
  dx[iter]=sort(x[j,(1:5)*3],decreasing = TRUE)[2]-sort(x[j,(1:5)*3])[2]
  iter = iter +1
}


this_fac = dy/(dx + 1.0)


return (this_fac)

}

##############################################

index_genes <- function(way){
  
  if (way==0){
    stds <- matrix(0,K,1);
    for (i in 1:(K)) {
      
      stds[i] <- sd(x[i,1:n])
      
    }
    
    i_stds <- order(stds, decreasing=TRUE)
    return (i_stds)
  }
  
  #################################
  if (way==1){
    changex <- matrix(0, K, 5)
    
    
    i <- 1
    for (j in seq(3,15, by=3)){
      
      thisx <- x[1:K,j]
      thisy <- abs(thisx-x[1:K,j-treat])
      
      changex[1:K,i] <- thisy/(thisx + 0.001)
      
      i= i+1
    }
    
    mean_changex <- matrix(0,K,1)
    
    for (i in 1:K) {
      
      mean_changex[i] <- mean(changex[i,1:5])
    }
    
    
    i_stds <- order(mean_changex, decreasing=TRUE)
    
    return (i_stds)
  }
  ####################################
  
  if (way==2){
    
    
    stds <- matrix(0,K,1);
    for (i in 1:(K)) {
      
      stds[i] <- sd(x[i,1:n])
      
    }
    
    i_stds <- order(stds, decreasing=TRUE)
    
    ###########################
    
    foo <- paste(substring(label,1, nchar(label)-1), 'ins', sep="")
    
    that1 <- paste(foo, "pos", treat_val, 'all', sep="_" )
    print(that1)
    
    
    assign(that1,common_gens(x,i_stds, K,treat,0, change_tol))
    
    this1 <- eval(parse(text = that1))
    #print(this1)
    
    ##########################
    
    that2 <- paste(foo, "neg", treat_val, 'all', sep="_" )
    print(that2)
    
    
    assign(that2,common_gens(x,i_stds, K,treat,1, change_tol))
    
    
    this2 <- eval(parse(text = that2))
    #print(this2)
    
    ########################
    
    
    that3 <- paste(foo, "null", treat_val, 'all', sep="_" )
    
    
    assign(that3,null_gens(x,i_stds, K,treat, change_tol))
    
    
    this3 <- eval(parse(text = that3))
    #print(this3)
    
    
    
    #########################
    
    
    i_stds0 <- c(this1,this2,this3)
    
    
    ##########
    
    changex <- matrix(0, K, 5)
    
    
    i <- 1
    for (j in seq(3,15, by=3)){
      
      thisx <- x[1:K,j]
      thisy <- abs(thisx-x[1:K,j-treat])
      
      changex[1:K,i] <- thisy/(thisx + 0.001)
      
      i= i+1
    }
    
    mean_changex <- matrix(0,K,1)
    
    for (i in 1:K) {
      
      mean_changex[i] <- mean(changex[i,1:5])
    }
    
    #########
    
    this_i_stds <- order(mean_changex[i_stds0], decreasing=TRUE)
    
    i_stds <- i_stds0[this_i_stds]
    
    
    return (i_stds)
    
  }
  
  #########################################
  
  if (way==3){
    changex <- matrix(0, K, 5)
    
    i <- 1
    for (j in seq(3,15, by=3)){
      
      thisx <- x[1:K,j]
      thisy <- abs(thisx-x[1:K,j-treat])
      
      changex[1:K,i] <- thisy
      
      i= i+1
    }
    
    sum_changex <- matrix(0,K,1)
    
    for (i in 1:K) {
      
      sum_changex[i] <- sum(changex[i,1:5])/(sort(x[i,(1:5)*3],decreasing = TRUE)[1]-sort(x[i,(1:5)*3])[1]+ 1.0)
      
    }
        
    i_stds <- order(sum_changex, decreasing=TRUE)
    
    out <- list(i_stds, sum_changex) 
    
    return (out)
    }
  
    
    if (way==4){
      
      ksdist <- matrix(0, K, 1)
      
      
      for (i in 1:K){

        thisks <- ks.test(x[i,(1:5)*3], x[i,(1:5)*3-treat])

        ksdist[i] <- unname(unlist(thisks[2]))

                
      }
      
      
      
      i_stds <- which(ksdist<0.01)
      
      out <- list(i_stds, ksdist[i_stds]) 
      
      return(out)
      
      }
    
    
   
  
  #####
  
  
}


##########################################




common_gens_glu <- function(x,i_stds, ngens, treat,par, tol) {
  
  # common gems for individuals wrt to a treatment
  
  i=1
  
  thisx <- x[i_stds[1:ngens],2]
  
  thisy <- thisx-x[i_stds[1:ngens],2-treat]
  #print(length(thisy[(-1)^par*thisy>0]))
  
  thisins <- which((-1)^par*thisy > tol*thisx)
  
  i= i+1
  
  for (j in seq(4,15, by=3)){
    
    thisx <- x[i_stds[1:ngens],j]
    thisy <- thisx-x[i_stds[1:ngens],j-treat]
    
    #print(length(thisy[(-1)^par*thisy>0]))
    
    newins <- which((-1)^par*thisy > tol*thisx)
    
    thisins <- thisins[which(thisins %in% newins)]
    
    i= i+1
  }
  return (i_stds[thisins])
}

################################################


myerr <- function(x,y,a,b,alpha){
  
  beta= 0.5*b/(alpha)
  
  thiserr <- max(abs(x-y))
  
  if (thiserr<b && max(x)<a){
    
    return (alpha)
    
  } else if (thiserr < 0.5*b && max(x)<beta){
    
    return (alpha)
    
  } else {
    
    return (thiserr/max(x))
    
  }
  
  
}

#################################################


latex_table <- function(thismat){
  
  digs =3
  
  for (i in 1:dim(thismat)[1]){
    
    thisline <- paste(toString(i),"&", substring(toString(thismat[i,1]),1,digs+2), sep=" ")
    
    for (j in 2:dim(thismat)[2]){
      
      thatline = paste("&", substring(toString(thismat[i,j]),1,digs+2), sep=" ")
      
      thisline <- paste(thisline, thatline, sep=" ")
      
    }
    
    print(thisline)
    
  }
  
}

######################################################

peaks <- function(thisx){
  
  data <- data.frame(x  = thisx)


  p <- ggplot(data, aes(x = x)) +  geom_density(alpha = 0.7) + scale_x_continuous(limits = c(0, 50))

  # Extract density data from the ggplot build
  density_data <- ggplot_build(p)$data[[1]]

  # Find the x value corresponding to the maximum density (the fitted curve peak)
  peak_x <- density_data$x[which.max(density_data$y)]
  
  return (peak_x)
  
}


################################################

group_gene_inds <- function(label, treat, G,dx,dy,i_stds,footol, this_fac){
  
  thistolx0 <- peaks(dx[i_stds]);
  thisF0 <- ecdf(dx[i_stds])
  dtolx0 <- -(unname(quantile(dx[i_stds],thisF0(thistolx0)+0.1))- thistolx0)/0.2*(footol-0.75)

  thistoly0 <- peaks(dy[i_stds]);
  thisF0 <- ecdf(dy[i_stds])
  dtoly0 <- -(unname(quantile(dy[i_stds],thisF0(thistoly0)+0.1))- thistoly0)/0.2*(footol-0.75)

  thistolr0 <- peaks(this_fac[i_stds]);
  thisF0 <- ecdf(this_fac[i_stds])
  dtolr0 <- -(unname(quantile(this_fac[i_stds],thisF0(thistolr0)+0.1))- thistolr0)/0.2*(footol-0.75)


  if (G=="G1") {
    
  ########## G1
    gene_inds <- intersect_gens(list(i_stds[which(dx[i_stds] < thistolx0 + dtolx0)], i_stds[which(dy[i_stds] < thistoly0 + dtoly0)]), 1:2)
  
  return (gene_inds)
  } 
    
  if (G=="G2") {
    
    ########### G2
  
  thistolx1 <- unname(quantile(dx[i_stds], footol))
  gene_inds <- intersect_gens(list(i_stds[which(dx[i_stds] > thistolx1)], i_stds[which(dy[i_stds] < thistoly0 + dtoly0 )]), 1:2)
  
  return (gene_inds)
  }
    
  if (G=="G4"){
    
    
   ########## G4
  
  thistoly1 <- unname(quantile(dy[i_stds], footol))
  gene_inds <- intersect_gens(list(i_stds[which(dx[i_stds] < thistolx0 + dtolx0)], i_stds[which(dy[i_stds] > thistoly1)]), 1:2)
  
   return (gene_inds)
    
  }

 
  if (G=="G3"){    
    
   ########## G3
  
  thistolx1 <- unname(quantile(dx[i_stds], footol))
  thistoly1 <- unname(quantile(dy[i_stds], footol))

  gene_inds <- intersect_gens(list(i_stds[which(dx[i_stds] > thistolx1)], i_stds[which(dy[i_stds] > thistoly1)], i_stds[which(this_fac[i_stds] < thistolr0 + dtolr0)]), 1:3)
  
   return (gene_inds)
    
  }


  
}

