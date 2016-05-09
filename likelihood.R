# COMS 4761: Project
# Students: Wenrui Huang; Vahe Galstyan

# INPUT
# K: number of states in the model
# M: the size of the alphabet (4 in the case of {A, G, T, C})
# bpSequence: the sequence of emitted nuleotides
# pi0: initial state cmf for the model
# A: the transition matrix for the model
# E: the emission matrix for the model

# OUTPUT
# Log[P(bpSequence | HMM)]: the log-probability of observing the emitted sequence,
#                           given the model parameters
# Log[P(bpSequence | HMM)} / length(bpSequence): the log-probability of observing a single
#                                                base-pair from the emitted sequence

# NOTES:
# log-transform is used throughout the calculations to avoid underflow issues
# special treatment is given to arithmetic operations with log-s (+, *, ...)

likelihood = function(K, M, bpSequence, pi0, A, E){  
  
  pi0.log <- log(pi0)
  A.log <- log(A)
  E.log <- log(E)
  SeqLength <- length(bpSequence)
  
  alpha.matrix.log <- matrix(rep(-Inf, SeqLength*K), nrow=K)    
  alpha.matrix.log[,1] <- E.log[, bpSequence[1]] + pi0.log
  
  for (t in 2:SeqLength){    
    for (i in 1:K){
      term1 <- E.log[i, bpSequence[t]]
      
      alphaPrime = alpha.matrix.log[, t-1] + A.log[i, ]
      
      alpha.matrix.log[i,t] = term1 + max(alphaPrime) + 
        log(sum(exp(alphaPrime - max(alphaPrime))))      
    }        
  }
  LogLikelihood <- max(alpha.matrix.log[,SeqLength]) + 
    log(sum(exp(alpha.matrix.log[,SeqLength] - max(alpha.matrix.log[,SeqLength]))))
  
  return( list(LogLikelihood, LogLikelihood/length(bpSequence)))
}