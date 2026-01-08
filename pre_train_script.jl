
########################################################### pre main loop
# activation function, training and testing individuals are selected

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

##########################################

@rget group label treat_val


##############################################################
#######################   loading DEGs
R"""

this <- paste(label,"limma_DEG_inds",treat_val, sep="_")

this <- paste(this,".csv",sep="")

df1 = read.table(this, header = TRUE, sep = ""); # csv file with the data
gene_inds <- df1[[1]]

ng <- length(gene_inds)

"""

@rget gene_inds 

@rget ng

########## or 

######### deserialize data

group_label=""
for g in group
    group_label= group_label*g
end

this_gene_inds = deserialize("gene_inds_"*label*treat_val*group_label)

gene_inds= this_gene_inds[1]

@rput gene_inds

ng = length(gene_inds)

########## or 

######### compute DEGs 

include("deg_extreme_groups.jl")
##############################################################

@rget gene_inds group_inds footol
@rget label treat treat_val group


#########################################################
##########################################################################
#########################################################

@rget N
N=Int(N)
ntrain = 4; 

allperts= digits.(0:(2^N-1), base=2, pad=N)
this= sum.(allperts)

perts=allperts[myfindall(x -> x == ntrain, this)]

baseinds=(1:5)*3


pertind=1
traininds= baseinds[Bool.(perts[pertind])]

testinds= baseinds[.!Bool.(perts[pertind])]

@rput traininds testinds

#######################################################################

R"""

dirvec <- c(0,1) ################################################## change dir here


traindx=data.frame(t(x[gene_inds,traininds-treat*dirvec[1]])) # default training temps

testdx=data.frame(t(x[gene_inds,testinds-treat*dirvec[1]]))

################################

traindy=data.frame(t(x[gene_inds,traininds-treat*dirvec[2]])) # default training temps

testdy=data.frame(t(x[gene_inds,testinds-treat*dirvec[2]]))


"""

@rget traindx testdx traindy testdy
@rget traininds testinds

trainX=copy(transpose(Matrix(traindx)[1:length(traininds),1:ng]));
trainY=copy(transpose(Matrix(traindy)[1:length(traininds),1:ng]));

testX=copy(transpose(Matrix(testdx)[1:length(testinds),1:ng]));
testY=copy(transpose(Matrix(testdy)[1:length(testinds),1:ng]));


thisY= copy(trainY);
thisX= copy(trainX);

##############################################

X=Interval.(zeros(ng,1));
Y=Interval.(zeros(ng,1));

for i=1:ng

    this= Interval.(minimum(thisX[i,:]), maximum(thisX[i,:]));
    X[i]= this[1]

end

for i=1:ng

    this= Interval(minimum(thisY[i,:]), maximum(thisY[i,:]))
    Y[i]= this[1]

end

#####################################################

coverX= Interval.(zeros(ng,ntrain))

for i =1:ng

    coverX[i,:]=gene_cover_1(trainX[i,:])'
    
end

coverY= Interval.(zeros(ng,ntrain))

for i =1:ng

    coverY[i,:]=gene_cover_1(trainY[i,:])'
    
end


#####################################################################################

R"""
library(neuralnet)
"""

tol0= 0.4;
@rput tol0;

global thisend=200
act_func!("linstep2")

input_size = size(trainX,1)
hidden_size = ng
output_size = 1
epochs = 10000

R"""
lrate <- 0.005
nsteps <- 1e6
"""
    
@rget lrate 
learning_rate = lrate


#####################################################################################
###################################  end of pre processing 
#####################################################################################