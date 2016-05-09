# COMS 4761: Project
# Students: Wenrui Huang; Vahe Galstyan

# INPUT
# K: number of states in the model
# M: the size of the alphabet (4 in the case of {A, G, T, C})
# bpSequence: the sequence of emitted nucleotides
# pi0: initial condition for the initial state cmf for the model
# A: initial condition for the transition matrix for the model
# E: initial condition for the emission matrix for the model
# eps: the minimum difference between the log-likelihoods of observing the bpSequence
#      in consecutive fwd/bwd iterations required for termination for the loop
# bwMax: the maximum number of fwd/bwd iterations in the EM algorithm
# Nlocal: number of local EM runs

# OUTPUT
# pi0.log: the log of the globally inferred initial state cmf of the model
# A.log: the log of the globally inferred transition matrix for the model
# E.log: the log of the globally inferred emission matrix for the model
# LogLikelihoods: the log likelihoods of oberserving the bpSequence at each fwd/bwd iteration

# NOTES:
# log-transform is used throughout the calculations to avoid underflow issues
# special treatment is given to arithmetic operations with log-s (+, *, ...)

source('bwLearning.R')
source('randomIC.R')
source('likelihood.R')

bwLearningGlobal = function(K, M, bpSequence, eps, bwMax, Nlocal){
  
  LogLikelihoods <- rep(-Inf, Nlocal)
  pi0.inf <- list()
  A.inf <- list()
  E.inf <- list()
  
  for (n in 1:Nlocal){
    
    # generate random initial condition
    IC <- randomIC(K, M)
    pi0 <- IC[[1]]
    A <- IC[[2]]
    E <- IC[[3]]
    
    # find a locally optimum solution    
    paramLocal <- bwLearning (K, M, pi0, A, E, bpSequence, eps, bwMax)
    
    LogLikelihoods[n] <- max(paramLocal[[4]])
    pi0.inf[[n]] <- exp(paramLocal[[1]])
    A.inf[[n]] <- exp(paramLocal[[2]])
    E.inf[[n]] <- exp(paramLocal[[3]])    
  }
  return ( list(pi0.inf, A.inf, E.inf, LogLikelihoods) )
}