# COMS 4761: Project
# Students: Wenrui Huang; Vahe Galstyan

# INPUT
# K: number of states in the model
# M: the size of the alphabet (4 in the case of {A, G, T, C})

# OUTPUT
# Uniformly random HMM parameters
# pi0: initial state cmf for the model
# A: the transition matrix for the model
# E: the emission matrix for the model

randomIC = function(K, M){
  
  # uniform randomization of the initial probability distribution of states
  pi0 <- c(runif(K, 0, 1))
  pi0 <- pi0/sum(pi0)
  
  # uniform randomization of the state transition matrix
  A <- matrix(runif(K^2, 0, 1), nrow=K, ncol=K)
  for (j in 1:K){
    A[,j] <- A[,j]/sum(A[,j])
  }
  
  E <- matrix(runif(K*M, 0, 1), nrow=K, ncol=M)
  for (i in 1:K){
    E[i,] <- E[i,]/sum(E[i,])
  }
  
  return (list(pi0, A, E))
}