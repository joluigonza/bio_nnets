########################################### deserialize data
@rget group label treat_val

############################################################

group_label=""
for g in group
    group_label= group_label*g
end

this_gene_inds = deserialize("gene_inds_"*label*treat_val*group_label)

##############################################################

group_label=""
for g in group
    group_label= group_label*g
end

######

part=""

part="limma"

forWs1 = deserialize("forWs1_"*label*treat_val*group_label*part)
forBs1 = deserialize("forBs1_"*label*treat_val*group_label*part)
forWs2 = deserialize("forWs2_"*label*treat_val*group_label*part)
forBs2 = deserialize("forBs2_"*label*treat_val*group_label*part)
forTols = deserialize("forTols_"*label*treat_val*group_label*part)
forEnds = deserialize("forEnds_"*label*treat_val*group_label*part)
forConv = deserialize("forConv_"*label*treat_val*group_label*part)

###################################################################################
####################################################################################
### if local parallization 

group_label=""
for g in group
    group_label= group_label*g
end

part="p0"

forWs1 = deserialize("forWs1_"*label*treat_val*group_label*part)
forBs1 = deserialize("forBs1_"*label*treat_val*group_label*part)
forWs2 = deserialize("forWs2_"*label*treat_val*group_label*part)
forBs2 = deserialize("forBs2_"*label*treat_val*group_label*part)
forTols = deserialize("forTols_"*label*treat_val*group_label*part)
forEnds = deserialize("forEnds_"*label*treat_val*group_label*part)
forConv = deserialize("forConv_"*label*treat_val*group_label*part)

ng= length(forBs2)
###################################

parts_num=3;
part0= Int(floor(ng/parts_num));

for i=1:(parts_num-1)
    mule_num=i;
    part="p"*string(mule_num)

    forWs1temp = deserialize("forWs1_"*label*treat_val*group_label*part)
    forBs1temp = deserialize("forBs1_"*label*treat_val*group_label*part)
    forWs2temp = deserialize("forWs2_"*label*treat_val*group_label*part)
    forBs2temp = deserialize("forBs2_"*label*treat_val*group_label*part)
    forTolstemp = deserialize("forTols_"*label*treat_val*group_label*part)
    forEndstemp = deserialize("forEnds_"*label*treat_val*group_label*part)
    forConvtemp = deserialize("forConv_"*label*treat_val*group_label*part)

    ###################################

    ind0= (1 + mule_num*part0)
    ind1= (mule_num + 1)*part0  + Int(floor((mule_num + 1)/parts_num))*(ng%parts_num)

    ############################################

    forWs1[ind0:ind1]= forWs1temp[ind0:ind1]
    
    forBs1[ind0:ind1]= forBs1temp[ind0:ind1]
    forWs2[ind0:ind1]= forWs2temp[ind0:ind1]
    forBs2[ind0:ind1]= forBs2temp[ind0:ind1]
    forTols[ind0:ind1]= forTolstemp[ind0:ind1]
    forEnds[ind0:ind1]= forEndstemp[ind0:ind1]
    forConv[ind0:ind1]= forConvtemp[ind0:ind1]

end

########################################################################