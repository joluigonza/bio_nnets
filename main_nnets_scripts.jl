
#########################


redoinds = myfindall(x -> x == 0.0, forConv)


#########################################################################
######################## Predictions and errors
#########################################################################

thisY= [trainY testY];
thisX= [trainX testX];

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

######################################################

coverX= Interval.(zeros(ng,N))

for i =1:ng

    
    coverX[i,:]=gene_cover_1([trainX[i,:]; testX[i]])'
    

end

coverY= Interval.(zeros(ng,N))

for i =1:ng
    
    coverY[i,:]=gene_cover_1([trainY[i,:]; testY[i]])'
    

end



########################################################################
R"""

thisgenes <- df[[1]][gene_inds[1:10]]

"""

@rget thisgenes

############################################################
################################# interval errors

REs=zeros(ng,1);
SEs=zeros(ng,1);

for k=1:ng

    trainin= copy(k)
    
    printstyled(k,color = :blue)
    println("")

    y=Y[trainin];

    thisy= myround_i(y,3);

    println("y: $thisy")

    ############################################
    
    
    z2= ∅;
    for i=1:ntrain
    a1, z1, a2, thisz2=forward_pass_i(coverX[:,i], forWs1[k], forBs1[k], forWs2[k], forBs2[k]);
    z2= z2[1] ∪ thisz2[1];
    end
    
        
    thisz2= myround_i(z2,3); 
    println("z2: $thisz2")  

    thissum1= int_dist(y,z2);
    println("RE:  $thissum1")

    REs[k]= thissum1

    o2= ∅;
    for i=1:ntrain
    ~, ~, ~, thiso2=forward_pass_i(coverY[:,i], forWs1[k], forBs1[k], forWs2[k], forBs2[k]);
    o2= o2[1] ∪ thiso2[1];
    end
      
        
    thiso2= myround_i(o2,3); 
    println("o2: $thiso2")  
    thissum2= int_dist(y,o2);
    println("SE:  $thissum2")

    SEs[k]= thissum2

    #println(thisgenes[k]*"&"*"$thisy"*"&"*"$thisz2"*"&"*"$thissum1"*"&"*"$thiso2"*"&"*"$thissum2")


    
end

thismean= round(mean(REs),digits=3)
thisvar= round(var(REs),digits=3)
println("($thismean,$thisvar)")


thismean= round(mean(SEs),digits=3)
thisvar= round(var(SEs),digits=3)
println("($thismean,$thisvar)")


#####################################################################################
####################################################### latex


for k=1:10

    trainin= copy(k)
    
    
    println("\\hline")
    
    println("")

    y=Y[trainin];

    thisy= myround_i(y,3);


    ############################################
    
    
    z2= ∅;
    for i=1:ntrain
    a1, z1, a2, thisz2=forward_pass_i(coverX[:,i], forWs1[k], forBs1[k], forWs2[k], forBs2[k]);
    z2= z2[1] ∪ thisz2[1];
    end
    
        
    thisz2= myround_i(z2,3); 
    

    thissum1= int_dist(y,z2);
    

    o2= ∅;
    for i=1:ntrain
    ~, ~, ~, thiso2=forward_pass_i(coverY[:,i], forWs1[k], forBs1[k], forWs2[k], forBs2[k]);
    o2= o2[1] ∪ thiso2[1];
    end
      
        
    thiso2= myround_i(o2,3); 
    

    thissum2= int_dist(y,o2);
    

    
    println(thisgenes[k]*"&"*"$thisy"*"&"*"$thisz2"*"&"*"$thissum1"*"&"*"$thiso2"*"&"*"$thissum2")


    println("")
    
    println("\\\\")

    
    println("")
    
end
    
############################################################
############### interpolation error 

for k=1:10

    trainin= copy(k)
    
    printstyled(k,color = :blue)
    println("")


    thisy= round(testY[k],digits=3);

    println("y: $thisy")

    ############################################
    
    
    a1, z1, a2, z2=forward_pass(testX, forWs1[k], forBs1[k], forWs2[k], forBs2[k]);
        
        
    thisz2= round(z2[1],digits=3); 
    println("z2: $thisz2")  

    thissum1= round(abs(thisy-thisz2)/(thisy),digits=3);

    
end
    

###############################################################################
############################### Interactions 
#####################################################################

using PyPlot
pygui(true)
#plot(x, y, color="red", linewidth=2.0, linestyle="--")


for k=1:ng
        
    thisabs = (abs.(forWs2[k])*abs.(forWs1[k]))
        
    plot(sort(thisabs[1,:]))
    
end

######################################## Top interactions

gen_minmax= zeros(ng,2,2);

for k=1:ng
    thisnorm = abs.(forWs2[k])*abs.(forWs1[k])

    gen_minmax[k,1,1]= minimum(thisnorm)
    gen_minmax[k,1,2]= argmin(thisnorm)[2]

    gen_minmax[k,2,1]= maximum(thisnorm)
    gen_minmax[k,2,2]= argmax(thisnorm)[2]


end


###############

topinters= (sortperm(gen_minmax[:,2,1],rev=true))[1:10]


Int.([topinters gen_minmax[topinters,2,2]])


@rget gene_inds

topgens = gene_inds[topinters]


topgens2 = gene_inds[Int.(gen_minmax[topinters,2,2])]

@rput topgens topgens2


#########################

R"""
    

    for (i in 1:10){
    print(i)
    print(paste(df[[1]][topgens[i]],"&", df[[1]][topgens2[i]], sep=" "))
    
    }


"""

############################################################
##################################### network entropy
##########################################################

thisng=size(forBs2,1);

totalplus=zeros(thisng,thisng);
totalminus=zeros(thisng,thisng); 

for k=1:thisng

    thisabs = transpose(abs.(forWs2[k])).*abs.(forWs1[k])
    this = transpose((forWs2[k])).*(forWs1[k])
    
    ###############

    thisnum=(abs.(this)).*(this .> 0)
    
    thistotal= sum(thisnum, dims=1) .+ 10^(-10)

    thisp= (thisnum) ./ (thistotal)

    thislogp = log.((abs.(this))).*(this .> 0) .- log.(thistotal).*(this .> 0)

    totalplus[k,:]= sum(thisp.*thislogp, dims=1)

    ##################

    thisnum=(abs.(this)).*(this .< 0)
    
    thistotal= sum(thisnum, dims=1) .+ 10^(-10)

    thisp= (thisnum) ./ (thistotal)

    thislogp = log.((abs.(this))).*(this .< 0) .- log.(thistotal).*(this .< 0)

    totalminus[k,:]= sum(thisp.*thislogp, dims=1)


end

##############################

totalentropy= sum(sum(-totalplus .- totalminus, dims=2))./(thisng^2*log(thisng-1))

round(totalentropy, digits=3)

#########################################
################## some tests 

myfindall2(x -> x>1, thisp)


myfindall2(x -> x>0, thislogp)

test= sum(thisp.*thislogp, dims=1)

myfindall2(x -> x>0, test)

println(test[:,1]')

findall( isnan.(totalminus) .==1)

test= thisp.*thislogp

findall( isnan.(test) .==1)

########################################## entropy of G4 on the rest 

thissumG4=sum(-totalplus[:,(group_inds[1]+group_inds[2]+1):end] , dims=2)[1:(group_inds[1]+group_inds[2])]

totalentropyG4= sum(thissumG4)./(thisng^2*log(thisng-1))

############################################################################
######################################### dynamical well-being v1 (outdated)
###########################################################################

Z2= Interval.(zeros(ng,1))
O2= Interval.(zeros(ng,1))
Sum1= zeros(ng,1);
Sum2= zeros(ng,1);
pertSum1= zeros(ng,1);
pertSum2= zeros(ng,1);

pertpara=0.1
nperts=20

for k=1:ng
    trainin= copy(k)
    
    printstyled(k,color = :blue)
    println("")

    y=Y[trainin];

    thisy= myround_i(y,3);

    println("y: $thisy")
    
    a1, z1, a2, z2=forward_pass_i(X, forWs1[k], forBs1[k], forWs2[k], forBs2[k]);

    Z2[k]=z2[1]
       
    ~, ~, ~, o2=forward_pass_i(Y, forWs1[k], forBs1[k], forWs2[k], forBs2[k]);

    O2[k]=o2[1]
    
    thisz2= myround_i(z2[1],3); 

    println("z2: $thisz2")  
    
    thissum1= round((abs(inf(y)-inf(z2[1])) + abs(sup(y)-sup(z2[1])))/(abs(sup(y))+abs(inf(y))) ,digits=3)  
    println(thissum1)    
    
    thiso2= myround_i(o2[1],3)
    println("o2: $thiso2")

    thissum2= round((abs(inf(y)-inf(o2[1])) + abs(sup(y)-sup(o2[1])))/(abs(sup(y))+abs(inf(y))), digits=3)
    println(thissum2)  

    Sum1[k]= thissum1
    Sum2[k]= thissum2
    
    for j=1:nperts

        pertWs1= forWs1[k] + pertpara*norm(forWs1[k],Inf).*rand(ng,ng)
        pertBs1= forBs1[k] + pertpara*norm(forBs1[k],Inf).*rand(ng,1)
        pertWs2= forWs2[k] + pertpara*norm(forWs2[k],Inf).*rand(1,ng)        
        pertBs2= forBs2[k] + pertpara*norm(forBs2[k],Inf).*rand(1,1)
            
        a1, z1, a2, z2=forward_pass_i(X, pertWs1,pertBs1,pertWs2, pertBs2);
        
        ~, ~, ~, o2=forward_pass_i(Y,  pertWs1,pertBs1,pertWs2, pertBs2);

        thissum1= int_dist(y,z2)
        thissum2= int_dist(y,o2)

        println(z2)
        println(thissum1)

        pertSum1[k]= pertSum1[k] + thissum1;
        
        pertSum2[k]= pertSum2[k] + thissum2;


    end

    
end

pertSum1=1/nperts*pertSum1

pertdist1 = sum((pertSum1 ./ Sum1)[1:group_inds[1]+group_inds[2]])/(group_inds[1]+group_inds[2])

println("pertdist1 : $pertdist1")

pertSum2=1/nperts*pertSum2

pertdist2 = sum((pertSum2 ./ Sum2)[1:group_inds[1]+group_inds[2]])/(group_inds[1]+group_inds[2])

println("pertdist2 : $pertdist2")


####################################################################

############# dynamical well-being v2

##############################################################3 


coreng= group_inds[1]+group_inds[2]


coreng= group_inds[1]


coreng= ng

########################################################

dymErrs= zeros(coreng,1)

pertpara=1.0*10^(-16)

for k=1:(coreng)

    println(k)
    trainin= copy(k)

    @rput trainin
    

    Ws1 = forWs1[k]
    Bs1 = forBs1[k]
    Ws2 = forWs2[k]
    Bs2 = forBs2[k]

    global thisend= forEnds[k]

    pertWs1= Ws1 + pertpara*norm(Ws1,Inf).*rand(ng,ng)
    pertBs1= Bs1 + pertpara*norm(Bs1,Inf).*rand(ng,1)
    pertWs2= Ws2 + pertpara*norm(Ws2,Inf).*rand(1,ng)        
    pertBs2= Bs2 + pertpara*norm(Bs2,Inf).*rand(1,1)


    #############################################

    trainy = trainY[k,:]'

    global nweights = ones(1, ntrain)
    global nweights[[argmin(trainy)[2], argmax(trainy)[2]]] .= 1.0


    for j=1:2

    pertWs1, pertBs1, pertWs2, pertBs2, conv, dloss, loss1= train_neural_network(trainX, trainy, hidden_size, output_size, epochs, learning_rate, pertWs1, pertBs1, pertWs2, pertBs2)
    
    end
    
    thiserr = norm(Ws1-pertWs1) + norm(Bs1-pertBs1) + norm(Ws2-pertWs2) + norm(Bs2-pertBs2)


    dymErrs[k]= thiserr

    
end


sum(dymErrs)/coreng


