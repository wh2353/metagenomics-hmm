# COMS 4761: Project
# Students: Wenrui Huang; Vahe Galstyan

source('bwLearningGlobal.R')
source('randomIC.R')
source('likelihood.R')
source('functions.R')

K <- 4
M <- 4

ref <- bp2int("data/ref") # a sample reference sequence from the Silva Database (AJ400267.1)
readMatch <- bp2int("data/match") # a read that was matched almost perfectly to the ref bacterium
                                  # using the Muther pipeline

readMismatch1 <- bp2int("data/read1") # read #1 that was assigned as a mismatch using the Muther pipeline
readMismatch2 <- bp2int("data/read2") # read #2 that was assigned as a mismatch using the Muther pipeline
readMismatch3 <- bp2int("data/read3") # read #3 that was assigned as a mismatch using the Muther pipeline

eps <- 10^(-3) # min difference in LogLikelihoods in the consequtive f/d iterations, 
               # required for the termination of the loop
bwMax <- 1000 # max number of EM iterations in the HMM learning algorithms
Nlocal <- 5 # number of local optimizations

# fwd/bwd learning of the HMM model for the reference sequence; Nlocal different times
paramLearned <- bwLearningGlobal (K, M, Ref, eps, bwMax, Nlocal)

# finding the local runs with the highest likelihood
localLogL <- paramLearned[[4]]
GlobalN <- which(localLogL == max(localLogL))[1]

pi0.inf <- paramLearned[[1]][[GlobalN]]
A.inf <- paramLearned[[2]][[GlobalN]]
E.inf <- paramLearned[[3]][[GlobalN]]

LMatch <- exp(likelihood(K, M, readMatch[[1]], pi0.inf, A.inf, E.inf)[[2]]) # 0.2624464
LMatch_prime <- exp(likelihood(K, M, readMatch[[1]], c(1,1,1,1)/4, A.inf, E.inf)[[2]]) # 0.2619738 (uniform pi0)

LMismatch1 <- exp(likelihood(K, M, readMismatch1[[1]], pi0.inf, A.inf, E.inf)[[2]]) # 0.2591945
LMismatch1_prime <- exp(likelihood(K, M, readMismatch1[[1]], c(1,1,1,1)/4, A.inf, E.inf)[[2]]) # 0.2587277 (uniform pi0)

LMismatch2 <- exp(likelihood(K, M, readMismatch2[[1]], pi0.inf, A.inf, E.inf)[[2]]) # 0.252411
LMismatch2_prime <- exp(likelihood(K, M, readMismatch2[[1]], c(1,1,1,1)/4, A.inf, E.inf)[[2]]) # 0.2519564 (uniform pi0)

LMismatch3 <- exp(likelihood(K, M, readMismatch3[[1]], pi0.inf, A.inf, E.inf)[[2]]) # 0.2598015
LMismatch3_prime <- exp(likelihood(K, M, readMismatch3[[1]], c(1,1,1,1)/4, A.inf, E.inf)[[2]]) # 0.2593336 (uniform pi0)