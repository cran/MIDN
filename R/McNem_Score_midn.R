McNem_Score_midn <-
function(alpha,SW,ppl,pmi,POWO)  {


# napprox_score

napprox_score <- function(alpha,ppl,pmi,POWO)  {
  
  sig_H <- sqrt(ppl + pmi)
  sig_A <- sqrt((ppl+pmi) - (ppl-pmi)**2)
  napprox <- (sig_H*qnorm(1-alpha) + sig_A*qnorm(POWO))**2 / (ppl-pmi)**2
  return(napprox)
}

# Ende napprox_score

#CritC_condNOC

CritC_condNOC <- function(N,ALPHA)  {
  
  PCeps <- 1/2
  CritC <- matrix(rep(0,(n+1)*2),ncol=2)
  CritC[1,1] <- 0
  CritC[1,2] <- ALPHA
  for (NOC in 1:N)
  { X <- NOC
  PX <- 1-pbinom(X,NOC,PCeps)
  while (PX <= ALPHA && X >= 1)
  { X <- X-1
  PX <- 1-pbinom(X,NOC,PCeps)  }
  C <- X+1
  PCPL <- 1-pbinom(C,NOC,PCeps)
  GAM <- (ALPHA-PCPL)/(PX-PCPL)
  CritC[NOC+1,1] <- C
  CritC[NOC+1,2] <- GAM      }
  return(CritC)
}

# Ende  CritC_condNOC

# POW_EX_NR

POW_EX_NR <- function(N,CritC,ppl,pmi)  {
  
  PBNOC <- c(rep(1,N+1))
  POWNR <- 0
  POC <- ppl+pmi
  P <- ppl/(ppl+pmi)
  PBNOC[1] <- pbinom(0,N,POC)
  for (NOC_ in 2:(N+1))
  { NOC <- NOC_-1
  PBNOC[NOC_] <- pbinom(NOC,N,POC) - pbinom(NOC-1,N,POC)    }
  for (NOC in 1:N)
  { C <- CritC[NOC+1,1]
  PBXGTC_NOC <- 1-pbinom(C,NOC,P)
  POWNR <- POWNR+PBXGTC_NOC*PBNOC[NOC+1]    }
  POW <- POWNR
  POW <- POW+CritC[1,2]*PBNOC[1]
  for (NOC in 1:N)
  { C <- CritC[NOC+1,1]
  gam <- CritC[NOC+1,2]
  PBXEQC_NOC <- pbinom(C,NOC,P) - pbinom(max(0,C-1),NOC,P)*sign(C)
  POW <- POW+gam*PBXEQC_NOC*PBNOC[NOC+1]     }
  POW_EX_NR <- matrix(c(0,0),ncol=2)
  POW_EX_NR[1,1] <- POW
  POW_EX_NR[1,2] <- POWNR
  return(POW_EX_NR)
}

# Ende  POW_EX_NR




nstart <- napprox_score(alpha,ppl,pmi,POWO)
n <- round(nstart)

CritC <- CritC_condNOC(n,alpha)
#CritC
pow_ex_nr <- POW_EX_NR(n,CritC,ppl,pmi)
POWEX <- pow_ex_nr[1,1]
nstart <- n

#pow_ex_nr


incr_n <- 1
incr5_n <- 5*incr_n

if (POWEX < POWO)
{ repeat
{ n <- n+incr5_n
CritC <- CritC_condNOC(n,alpha)
pow_ex_nr <- POW_EX_NR(n,CritC,ppl,pmi)
POWEX <- pow_ex_nr[1,1]
if (POWEX >= POWO) break                }
  
  n <- n-incr5_n
  
  repeat
  { n <- n+incr_n
  CritC <- CritC_condNOC(n,alpha)
  pow_ex_nr <- POW_EX_NR(n,CritC,ppl,pmi)
  POWEX <- pow_ex_nr[1,1]
  if (POWEX >= POWO) break                }   }  else
  { repeat
  { n <- n-incr5_n
  CritC <- CritC_condNOC(n,alpha)
  pow_ex_nr <- POW_EX_NR(n,CritC,ppl,pmi)
  POWEX <- pow_ex_nr[1,1]
  if (POWEX < POWO) break                }
    repeat
    { n <- n+incr_n
    CritC <- CritC_condNOC(n,alpha)
    pow_ex_nr <- POW_EX_NR(n,CritC,ppl,pmi)
    POWEX <- pow_ex_nr[1,1]
    if (POWEX >= POWO) break                }   }

Nex <- n

POWNR <- pow_ex_nr[1,2]

if (POWNR < POWO)
{ repeat
{ n <- n+incr5_n
CritC <- CritC_condNOC(n,alpha)
pow_ex_nr <- POW_EX_NR(n,CritC,ppl,pmi)
POWNR <- pow_ex_nr[1,2]
if (POWNR >= POWO) break                 }
  n <- n-incr5_n
  repeat
  { n <- n+incr_n
  CritC <- CritC_condNOC(n,alpha)
  pow_ex_nr <- POW_EX_NR(n,CritC,ppl,pmi)
  POWNR <- pow_ex_nr[1,2]
  if (POWNR >= POWO) break                 }  }  else
  { repeat
  { n <- n-incr5_n
  CritC <- CritC_condNOC(n,alpha)
  pow_ex_nr <- POW_EX_NR(n,CritC,ppl,pmi)
  POWNR <- pow_ex_nr[1,2]
  if (POWNR < POWO) break                 }
    repeat
    { n <- n+incr_n
    CritC <- CritC_condNOC(n,alpha)
    pow_ex_nr <- POW_EX_NR(n,CritC,ppl,pmi)
    POWNR <- pow_ex_nr[1,2]
    if (POWNR >= POWO) break                 }  }

Nnr <- n

mid_n <- round((Nex+Nnr)/2 + 10**-6)


erg <- c(nstart,Nex,POWEX,Nnr,POWNR,mid_n)

}
