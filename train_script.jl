
#####################################################################################
########################################################### main loop
############################################################################################


forWs1 = Array{Any}(undef,ng);
forBs1 = Array{Any}(undef,ng);
forWs2 = Array{Any}(undef,ng);
forBs2 = Array{Any}(undef,ng);

forTols = zeros(ng,3);
forEnds = zeros(ng,1);
forConv = zeros(ng,1);
    

######################################################
### local parallelization

# mule_num=0;

# parts_num=3;

# part0= Int(floor(ng/parts_num))

# ind0= (1 + mule_num*part0)
# ind1= (mule_num + 1)*part0  + Int(floor((mule_num + 1)/parts_num))*(ng%parts_num)

###############################################
ind0=0;
ind1=ng;

for trainin= ind0: ind1
    
    #global thisend=200*ones(ng)  # restart thisend here in case we want to optimze it for each network
    
    global thisend=200  # restart thisend here in case we want to optimze it for each network
    
    global nweights= ones(1,ntrain)

    println("")
    println("trainin: $trainin")

    try
               
    #########################   

    @rput trainin
    
    ########################
    
    R"""       
   
    thistrainy=as.matrix(unname(traindy[trainin]))
    
    thistesty=as.matrix(unname(testdy[trainin]))


    thiserr <- "1.0"
    
    counter <- 0

    while ( !is.numeric(thiserr) && counter < 100) {
    
      
    try({
    nn <- neuralnet(thistrainy ~., data.frame(traindx, thistrainy), hidden=c(ng), act.fct="logistic",  linear.output = TRUE, algorithm="backprop", learningrate = lrate, stepmax = nsteps)    
    
    thiserr <-  max(abs(thistesty-predict(nn, testdx)))/(max(abs(thistesty))+10^(-10))
       
    }, silent=TRUE)

    counter <- counter +1
    print(paste("counter: ", as.character(counter)))
    print(paste("error: ", as.character(thiserr)))
    
    }


    """
##########################
@rget thiserr


if !(thiserr isa Number)  
    
    println("trainin: $trainin")
    println("error. no initial convergence")
    break  
    
end

tol0 = copy(thiserr);
@rput tol0

##########################

    R"""
    ws1 <- nn$weights[[1]][1]
    ws1 <- ws1[[1]]
    
    ws2 <- nn$weights[[1]][2]
    ws2 <- ws2[[1]]

    predy <- predict(nn, testdx)
    """

    @rget thistrainy thistesty predy 

    @rget ws1 ws2


    @rget trainin dirvec label
    trainin = Int(trainin)

    ##############################################
    
    @rget traininds testinds

    trainy=copy(transpose(thistrainy));

    testy=copy(transpose(thistesty));

    ##################################################


    W1_0= ws1[2:end,:];
    
    b1_0= ws1[1,:];

    W2_0= transpose(ws2[2:end]);
    b2_0= ws2[1];

    #################################
      
    global nweights[[argmin(trainy)[2], argmax(trainy)[2]]] .= 1.0

    #####################################

    W1, b1, W2, b2, conv, dloss, loss1, loss2= train_neural_network_0(trainX, trainy, hidden_size, output_size, epochs, learning_rate, transpose(W1_0), b1_0, W2_0, b2_0);

    
    #####################################

    for k=1:4
        println("continuation: $k")
        W1, b1, W2, b2,conv,~,~= train_neural_network(trainX, trainy, hidden_size, output_size, epochs, learning_rate, W1, b1, W2, b2);
    end

    #################################################
    thisWs1=Array{Any}(undef,1);
    thisBs1=Array{Any}(undef,1);
    thisWs2=Array{Any}(undef,1);
    thisBs2=Array{Any}(undef,1);
    

    thisWs1[1]= W1;
    thisBs1[1]= b1;
    thisWs2[1]= W2;
    thisBs2[1]= b2;

    ##############################################
    ################################################################### improving loop
    
    iters = 10; # number or net candidates
    @rput iters

    thisTols=zeros(iters+1,1);
    
    thisEnds=zeros(1,iters+1);
    thisEnds[1]=thisend;

    thisTols[1]=thiserr;
    
    thisConv=zeros(iters+1,1);
    thisConv[1]=conv;

    for j=1:iters

        println("")
        println("trainin: $trainin")
        println("iters: $j")
            
                
        # global thisend=200*ones(ng)  # restart thisend here in case we want to optimze it for each network

        global thisend =200 # restart thisend here in case we want to optimze it for each network
        
        global nweights= ones(1,ntrain)

    
        R"""       
    
        thiserr <- "1.0"
        
        counter <- 0

        while ( !is.numeric(thiserr) && counter < 100) {
        
        
        try({
        nn <- neuralnet(thistrainy ~., data.frame(traindx, thistrainy), hidden=c(ng), act.fct="logistic",  linear.output = TRUE, algorithm="backprop", learningrate = lrate, stepmax = nsteps)    
        
        thiserr <-  max(abs(thistesty-predict(nn, testdx)))/(max(abs(thistesty))+10^(-10))
        
            
        }, silent=TRUE)

        counter <- counter +1
        print(paste("counter: ", as.character(counter)))
        print(paste("error: ", as.character(thiserr)))
        
        }


        """
        ##########################

        @rget thiserr


        if !(thiserr isa Number)  # counter is 100
            
            
            println("error. no initial convergence")
            println("trainin: $trainin")

            break  
            
        end

        tol0 = copy(thiserr);
        @rput tol0

        ##########################

        R"""
        ws1 <- nn$weights[[1]][1]
        ws1 <- ws1[[1]]
        
        ws2 <- nn$weights[[1]][2]
        ws2 <- ws2[[1]]

        predy <- predict(nn, testdx)
        """

        @rget predy 

        @rget ws1 ws2


        ##############################################
        
        W1_0= ws1[2:end,:];
        
        b1_0= ws1[1,:];

        W2_0= transpose(ws2[2:end]);
        b2_0= ws2[1];

        #################################
            
        global nweights[[argmin(trainy)[2], argmax(trainy)[2]]] .= 1.0

        #####################################

        W1, b1, W2, b2, conv, dloss, loss1, loss2= train_neural_network_0(trainX, trainy, hidden_size, output_size, epochs, learning_rate, transpose(W1_0), b1_0, W2_0, b2_0);

        #################################################
        
        for k=1:4
            println("continuation: $k")
            W1, b1, W2, b2,conv,~,~= train_neural_network(trainX, trainy, hidden_size, output_size, epochs, learning_rate, W1, b1, W2, b2);
        end


        ###############################################

        push!(thisWs1,W1)
        push!(thisBs1,b1)
        push!(thisWs2,W2)
        push!(thisBs2,b2)

        thisTols[j+1]=thiserr;
        #thisEnds[:,j+1]=thisend;
        thisEnds[j+1]=thisend;
        thisConv[j+1]=conv;


    end ########################################## end of trainings loop (improving loops)

    println("trainin: $trainin")

    #####################################################
    ################################################################################
    #######################################################################

    y=Y[trainin];

    thisIns=zeros(size(thisBs2,1),1);
    
    thisOut=zeros(size(thisBs2,1),1);


    for k=1:size(thisBs2,1)
        println(k)

        z2= ∅;
        for i=1:ntrain
        a1, z1, a2, thisz2=forward_pass_i(coverX[:,i], thisWs1[k], thisBs1[k], thisWs2[k], thisBs2[k]);
        z2= z2[1] ∪ thisz2[1];
        end
        
        println("z2: $(z2[1])")        
         
        thisIns[k]=int_dist(y,z2);

        o2= ∅;
        for i=1:ntrain
        ~, ~, ~, thiso2=forward_pass_i(coverY[:,i], thisWs1[k], thisBs1[k], thisWs2[k], thisBs2[k]);
        o2= o2[1] ∪ thiso2[1];
        end
        
        println("o2: $(o2[1])")
       
        thisOut[k]=int_dist(y,o2);

    end


    thisarg= argmin(thisIns)[1]
    println("")
    printstyled(thisarg; color = :green)
    println("")
    printstyled(thisIns[thisarg]; color = :green)
    println("")
    
    printstyled(y; color = :blue)

    #################################################

    k=argmin(thisIns)[1]

    forWs1[trainin]=thisWs1[k]
    forBs1[trainin]=thisBs1[k]
    forWs2[trainin]=thisWs2[k]
    forBs2[trainin]=thisBs2[k]

    forTols[trainin,1:3]=[thisTols[k], thisIns[k], thisOut[k]]
    forEnds[trainin]= thisEnds[k]
    forConv[trainin]= thisConv[k]


    catch
    println("error somewhere")
    println("trainin: $trainin")
    break

    end
        
        

end

################################################
####################################################################################################
##################################################### end main loop