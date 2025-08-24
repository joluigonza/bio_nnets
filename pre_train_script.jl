
########################################################### pre main loop
# activation function, training and testing individuals are selected

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


act_func!("linstep2")

input_size = size(trainX,1)
hidden_size = ng
output_size = size(trainy,1)
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