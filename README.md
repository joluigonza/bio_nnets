# bio_nnets

1) loading networks from file

run pkgs_files_setup.jl  (load pkgs and function files)
run pre_train_script.jl (load data and  genes, activation function selected)
run loading_scripts.jl  (load networks)


2) training networks:

setup the nnet_setup.jl file

run pkgs_files_setup.jl   (load pkgs and function files)

run pre_train_script.jl (load data and  genes, activation function selected)
if needed run deg_extreme_group.jl (compute new DEGs)

run train_script.jl (train network) 
