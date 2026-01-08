
R"""

df = read.table("TPM_5-species.csv", header = TRUE, sep = ""); # csv file with the data

K <- nrow(df) #number of transcripts 
n <- 15 # 5 individuals * 3 temperatures
N <- 5 # individuals

"""

###########################################

R"""

# xloc <- (1:15)+1 #camel
# label <- 'camel'

xloc <- c(19,20,21,24,25,26,29,30,31,34,35,36,39,40,41) #human
label <- 'humans'

# xloc <- 57:71 # rhino
# label <- 'rhinos'

#############################################

treat <- 1  # 1 is 41, 2 is 32
treat_val <- "41"

# treat <- 2  # 1 is 41, 2 is 32
# treat_val <- "32"

### change_tol <- 0.1


#group <- c("G1")


group <- c("G1","G2","G3")


#group <- c("G1","G2","G4")

#group <- c("G3")


"""
##########################################
########################################################################

# R"""
# df = read.table("TPM_9-species.tsv", header = TRUE, sep = "");

# K <- nrow(df) #number of transcripts 
# n <- 15 # 5 individuals * 3 temperatures
# N <- 5 # individuals

# # S is 18 
# xloc <- c(18,19,20,23,24,25,28,29,30,33,34,35,38,39,40) #tpm_9 file
# label <- 'bats'


# treat <- 1  # 1 is 41, 2 is 32
# treat_val <- "41"

# # treat <- 2  # 1 is 41, 2 is 32
# # treat_val <- "32"

# ### change_tol <- 0.1

# group <- c("G1","G2","G4")

# """

