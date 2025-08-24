
group_label=""
for g in group
    group_label= group_label*g
end


this_gene_inds= Array{Any}(undef,3)
this_gene_inds[1]= gene_inds
this_gene_inds[2]= group_inds
this_gene_inds[3]= footol

serialize("gene_inds_"*label*treat_val*group_label, this_gene_inds)

############################# final save

group_label=""
for g in group
    group_label= group_label*g
end


part=""

serialize("forWs1_"*label*treat_val*group_label*part, forWs1)
serialize("forBs1_"*label*treat_val*group_label*part, forBs1)
serialize("forWs2_"*label*treat_val*group_label*part, forWs2)
serialize("forBs2_"*label*treat_val*group_label*part, forBs2)
serialize("forTols_"*label*treat_val*group_label*part, forTols)
serialize("forEnds_"*label*treat_val*group_label*part, forEnds)
serialize("forConv_"*label*treat_val*group_label*part, forConv)

####################################
### if local paralellization


group_label=""
for g in group
    group_label= group_label*g
end


parts_num=3;
part0= Int(floor(ng/parts_num));


part="p"*string(mule_num)

serialize("forWs1_"*label*treat_val*group_label*part, forWs1)
serialize("forBs1_"*label*treat_val*group_label*part, forBs1)
serialize("forWs2_"*label*treat_val*group_label*part, forWs2)
serialize("forBs2_"*label*treat_val*group_label*part, forBs2)
serialize("forTols_"*label*treat_val*group_label*part, forTols)
serialize("forEnds_"*label*treat_val*group_label*part, forEnds)
serialize("forConv_"*label*treat_val*group_label*part, forConv)

#####################################
### saving to csv

##################### saving


R"""

this <- data.frame(df[[1]][gene_inds])

write.table(this, "humans_g1g2g4_32.csv", col.names=FALSE, row.names=FALSE)

#write.table(this, "humans_G3_32_v1.csv", col.names=FALSE, row.names=FALSE)


"""


R"""

this <- data.frame(df[[1]][i_stds])

write.table(this, "humans_DEG_41.csv", col.names=FALSE, row.names=FALSE)


"""

#####################
