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
#      in consecutive fwd/bwd iterations required for termination of the loop
# bwMax: the maximum number of fwd/bwd iterations in the EM algorithm

# OUTPUT
# pi0.log: the log of the inferred initial state cmf of the model
# A.log: the log of the inferred transition matrix for the model
# E.log: the log of the inferred emission matrix for the model
# LogLikelihoods: the log likelihoods of oberserving the bpSequence at each fwd/bwd iteration

# NOTES:
# log-transform is used throughout the calculations to avoid underflow issues
# special treatment is given to arithmetic operations with log-s (+, *, ...)

bwLearning = function(K, M, pi0, A, E, bpSequence, eps, bwMax){      
  
  pi0.log <- log(pi0)
  A.log <- log(A)
  E.log <- log(E)
  SeqLength <- length(bpSequence)      
  
  # ------- bw iterations ------
  LogLikelihoods <- rep(0, times=bwMax)
    
  for (bw in 1:bwMax){      

    # -------------------- Log[alpha] coefficients ------------------------
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
    LogLikelihoods[bw] <- max(alpha.matrix.log[,SeqLength]) + 
      log(sum(exp(alpha.matrix.log[,SeqLength] - max(alpha.matrix.log[,SeqLength]))))
            
    # -------------------- Log[beta] coefficients ------------------------  
    beta.matrix.log <- matrix(rep(-Inf, SeqLength*K), nrow=K)
    beta.matrix.log[,SeqLength] = 0
    
    for (t in (SeqLength-1):1){
      for (i in 1:K){                      
        terms = beta.matrix.log[, t+1] + E.log[, bpSequence[t+1]] + A.log[, i]
        if (max(terms) == -Inf){
          beta.matrix.log[i,t] = -Inf
        }
        if (max(terms) != -Inf){
          beta.matrix.log[i,t] = max(terms) + log(sum(exp(terms - max(terms))))          
        }
      }
    }
            
    # -------------------- Log[P(S(t) | bpSequence)] ------------------------
    PSt.log <- alpha.matrix.log + beta.matrix.log - LogLikelihoods[bw]
    
    # -------------------- Log[P(S(t), S(t-1) | bpSequence)] ------------------------
    PSSt.log <- array(rep(-Inf, K^(2)), dim = c(K, K, SeqLength))
    for (t in 2:SeqLength){
      for (j in 1:K){       
        for (l in 1:K){
          PSSt.log[l, j, t] <- alpha.matrix.log[j, t-1] +
            E.log[l,bpSequence[t]] + A.log[l, j] + beta.matrix.log[l, t] - LogLikelihoods[bw]
        }
      }
    }
    
    
    ##################### Maximization Steps #####################
    
    # -------------------- Maximization: pi0 ------------------------  
    for (i in 1:K){
      terms <- PSt.log[,1]  
      pi0.log[i] <- PSt.log[i,1] - (max(terms) + log(sum(exp(terms - max(terms)))))
    }
    
    # -------------------- Maximization: A ------------------------  
    for (m in 1:K){
      for (n in 1:K){        
        terms <- PSSt.log[m,n,2:SeqLength]
        
        if (max(terms) == -Inf){
          A.log[m,n] <- -Inf
        }
        
        if (max(terms) != -Inf){
          A.log[m,n] <- max(terms) + log(sum(exp(terms - max(terms))))
        }                
      }
    }
    for (n in 1:K){
      A.log[,n] <- A.log[,n] - max(A.log[,n]) - log(sum(exp(A.log[,n] - max(A.log[,n]))))
    }
    
    # -------------------- Maximization: E ------------------------  
    
    for (k in 1:K){
      for (m in 1:M){
        terms <- PSt.log[k, which(bpSequence == m)]
        
        if (max(terms) == -Inf){
          E.log[k, m] <- -Inf
        }
        
        if (max(terms) != -Inf){
          E.log[k, m] <- max(terms) + log(sum(exp(terms - max(terms))))  
        }
      }
    }
    
    for (k in 1:K){
      E.log[k,] <- E.log[k,] - max(E.log[k,]) - log(sum(exp(E.log[k,] - max(E.log[k,]))))
    }     

    # -------------------- Convergence ------------------------    
#     if (bw > 1){
#       print(LogLikelihoods[bw]-LogLikelihoods[bw-1])      
#     }               
    
    if (bw > 1 && LogLikelihoods[bw]-LogLikelihoods[bw-1] <= eps){
      LogLikelihoods <- LogLikelihoods[1:bw]
      break
    }           
    
  }
  
  return (list(pi0.log, A.log, E.log, LogLikelihoods))
}