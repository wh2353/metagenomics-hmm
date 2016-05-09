# COMS 4761: Project
# Students: Wenrui Huang; Vahe Galstyan

# Used to check the consistency of the bwLearning algorithm
# The condition of consistency is that the likelihood of observing the synthetic
# nucleotide sequence under the inferred parameters should be higher than that under the
# initially chosen parameters.

source('bwLearningGlobal.R')
source('randomIC.R')
source('likelihood.R')

K <- 4
M <- 4

IC <- randomIC(K, M)
pi0 <- IC[[1]]
A <- IC[[2]]
E <- IC[[3]]

Nbp <- 200
stateSequence <- rep(0, Nbp)
bpSequence <- rep(0, Nbp)

stateSequence[1] <- sample(1:K, size=1, prob=pi0)

for (t in 2:Nbp){
  stateSequence[t] <- sample(1:K, size=1, prob = A[,stateSequence[t-1]])
}

for (t in 1:Nbp){
  bpSequence[t] <- sample(1:M, size=1, prob = E[stateSequence[t],])
}

eps <- 10^(-3)
bwMax <- 1000
Nlocal <- 2

paramLearned <- bwLearningGlobal (K, M, bpSequence, eps, bwMax, Nlocal)

localLogL <- paramLearned[[4]]
GlobalN <- which(localLogL == max(localLogL))[1]

pi0.inf <- paramLearned[[1]][[GlobalN]]
A.inf <- paramLearned[[2]][[GlobalN]]
E.inf <- paramLearned[[3]][[GlobalN]]

LogL0 <- likelihood(K, M, bpSequence, pi0, A, E)[[1]]
LogL.inf <- likelihood(K, M, bpSequence, pi0.inf, A.inf, E.inf)[[1]]

if (LogL.inf < LogL0){
  print("FAIL")
}
else{
  print("PASS")
}