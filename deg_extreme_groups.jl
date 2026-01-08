
include("nnet_setup.jl")

##############################################

R"""

x=matrix(0,K,n);

# find the location of the data in nnet_setup.jl
for (i in 1:(K)) {
  thisin=1
  for (j in xloc) {
    x[i,(thisin)]=df[[j]][i]
    thisin=thisin+1
  }
}
"""



R"""

thistol <- (N-2) # number of individuals - some tolerance


changex <- matrix(0, K, n)

i <- 1
for (j in seq(3,15, by=3)){
  
  thisx <- x[1:K,j]
  thisy <- abs(thisx-x[1:K,j-treat])
  
  changex[1:K,i] <- thisy
  
  i= i+1
}

sum_changex <- matrix(0,K,1)
dy=matrix(0,K,1);
dx=matrix(0,K,1);

for (i in 1:K) {

    thisdx= sort(x[i,(1:5)*3],decreasing = TRUE)[2]-sort(x[i,(1:5)*3])[2]
    thisdy= sort(x[i,(1:5)*3-treat],decreasing = TRUE)[2]-sort(x[i,(1:5)*3-treat])[2]    


    dx[i]=thisdx
    dy[i]=thisdy
    sum_changex[i] <- sum(sort(changex[i,1:5])[2:4])/(thisdx + thistol/thistol)*(((sum(x[i,(1:5)*3]>0)>1)+(sum(x[i,(1:5)*3-treat]>0)>1))>0)*1.0
    
  
  
}


i_stds <- order(sum_changex, decreasing=TRUE)


this <- sum_changex

thisng <- length(this[this>(1.5*thistol)])
print(thisng)


i_stds <- i_stds[1:thisng]


this_fac <- (dy)/(dx+ 10^(-10)); # avoid dividing by 0

"""

############################################################


R"""
footol <- 0.88 #0.85   #0.75

gene_inds <- vector()
group_inds <- vector()

for (thisg in group){

    thisgene_inds <- group_gene_inds(label,treat,thisg,dx,dy,i_stds,footol,this_fac)
    group_inds <- c(group_inds, length(thisgene_inds))
    gene_inds <- c(gene_inds, thisgene_inds)
    
} 


ng <- length(gene_inds)


"""


@rget  ng


###############################
