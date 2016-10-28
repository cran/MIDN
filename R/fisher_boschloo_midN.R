fisher_boschloo_midN <-
function(alpha,SW,p1,p2,POWO,mton_a,mton_b)  {

mton_ab <- mton_a + mton_b
lambda <- mton_a / mton_ab

# Napprox_score

Napprox_score <- function(alpha,p1,p2,POWO,lambda) {
  
  
  p0 <- lambda*p1 + (1-lambda)*p2
  sig_H <- sqrt(p0*(1-p0)*(1/lambda + 1/(1-lambda)))
  sig_A <- sqrt(p1*(1-p1)/lambda + p2*(1-p2)/(1-lambda))
  napprox <- (sig_H*qnorm(1-alpha) + sig_A*qnorm(POWO))**2/(p1-p2)**2
  napprox <- ceiling(napprox)
  N1 <- round(lambda*napprox + 10**-6)
  N2 <- round((1-lambda)/lambda*N1 + 10**-6)
  Napprox <- matrix(c(0,0),ncol=2) 
  Napprox[1,1] <- N1
  Napprox[1,2] <- N2
  return(Napprox)
}  

# Ende Napprox_score

#CritC_condS

CritC_condS <- function(M,N,alpha)  {
  
  RHO_0 <- 1
  NN <- M + N
  
  CritC <- matrix(rep(0,(NN*2+2)),ncol=2)
  CritC[1,1] <- 0
  CritC[1,2] <- alpha
  for (S in 1:(NN-1))
  { X <- min(S,M)
  PX <- 1-pFNCHypergeo(X-1,M,N,S,RHO_0)
  while (PX <= alpha && X >= max(0,S-N)+1)
  { X <- X-1
  PX <- 1
  if (X > max(0,S-N))
    PX <- 1-pFNCHypergeo(X-1,M,N,S,RHO_0)  }
  C <- X 
  PCPL <- 1-pFNCHypergeo(C,M,N,S,RHO_0)
  GAM <- (alpha-PCPL)/(PX-PCPL)
  CritC[S+1,1] <- C
  CritC[S+1,2] <- GAM    }
  CritC[NN+1,1] <- M
  CritC[NN+1,2] <- alpha
  return(CritC)
}

# Ende CritC_condS

# CritfctXY

CritfctXY <- function(M,N,CritC)   {
  
  phi <- matrix(rep(0,((M+1)*(N+1))),ncol=N+1)
  for (x in 0:M)
  { for (y in 0:N)
  { s <- x+y
  if (x < CritC[s+1,1])
    phi[x+1,y+1] <- 0
  if (x == CritC[s+1,1])
    phi[x+1,y+1] <- CritC[s+1,2]
  if (x > CritC[s+1,1])
    phi[x+1,y+1] <- 1     }    }
  return(phi) 
}

# Ende CritfctXY

# KX_Y

KX_Y <- function(M,N,phi)   {
  
  KX_Y <- c(rep(1,N+1))
  phi <- floor(phi)
  for (Y_ in 1:(N+1))
  { X_ <- 1
  while (phi[X_,Y_] == 0 && X_ <= M)
  { X_ <- X_+1  }
  KX_Y[Y_] <- X_-1
  if (KX_Y[Y_] == M && phi[M+1,Y_] == 0)
  { KX_Y[Y_] <- M+1  }  }
  return(KX_Y)
}

# Ende KX_Y 

# POW_EX_NR

POW_EX_NR <- function(M,N,KX_Y,CritC,P1,P2)   {
  
  PBY <- c(rep(1,N+1))
  POWNR <- 0
  pi2 <- P2
  PBY[1] <- pbinom(0,N,pi2)
  for (Y_ in 2:(N+1))
  {  Y <- Y_-1
  PBY[Y_] <- pbinom(Y,N,pi2) - pbinom(Y-1,N,pi2)   }
  
  pi1 <- P1
  for (Y_ in 1:(N+1))
  { Y <- Y_-1
  if (KX_Y[Y_] == 0) PBXGEK_Y <- 1
  if (KX_Y[Y_] == M+1) PBXGEK_Y <- 0
  if (KX_Y[Y_] >= 1 && KX_Y[Y_] <= M)
  { PBXGEK_Y <- 1-pbinom(KX_Y[Y_]-1,M,pi1)  }
  
  PBYEQY <- PBY[Y_]
  if (min(PBXGEK_Y,PBYEQY) <= 0) INCR <- 0
  if (min(PBXGEK_Y,PBYEQY) > 0)
  { LPBX <- log(PBXGEK_Y)
  LPBY <- log(PBYEQY)
  INCR <- exp(LPBX+LPBY)  }
  POWNR <- POWNR + INCR          }
  
  POW <- POWNR
  POW <- POW + CritC[1,2]*pbinom(0,M,pi1)*pbinom(0,N,pi2)
  NN <- M+N
  for (S in 1:NN)
  { x0 <- CritC[S+1,1]
  gam <- CritC[S+1,2]
  y0 <- S-x0
  POW <- POW + gam*(pbinom(x0,M,pi1)-pbinom(max(0,x0-1),M,pi1)*sign(x0)) *
    (pbinom(y0,N,pi2)-pbinom(max(0,y0-1),N,pi2)*sign(y0))  }
  POW_EX_NR <- matrix(c(0,0),ncol=2)
  POW_EX_NR[1,1] <- POW
  POW_EX_NR[1,2] <- POWNR
  return(POW_EX_NR)
}

# Ende POW_EX_NR









Nstart <- Napprox_score(alpha,p1,p2,POWO,lambda)
Nstart

m <- Nstart[1,1]
n <- Nstart[1,2]

m <- round((m+n)/mton_ab + 10**-6) * mton_a
n <- round((m+n)/mton_ab + 10**-6) * mton_b

CritC <- CritC_condS(m,n,alpha)

phi <- CritfctXY(m,n,CritC)

kx_y <- KX_Y(m,n,phi)

pow_ex_nr <- POW_EX_NR(m,n,kx_y,CritC,p1,p2)
POWEX <- pow_ex_nr[1,1]

mstart <- m
nstart <- n
s_mstart <- mstart
s_nstart <- nstart


incr_m <- mton_a
incr_n <- mton_b
incr5_m <- 5*incr_m
incr5_n <- 5*incr_n

if (POWEX < POWO)
{ repeat
{ m <- m+incr5_m
n <- n+incr5_n
CritC <- CritC_condS(m,n,alpha)
phi <- CritfctXY(m,n,CritC)
kx_y <- KX_Y(m,n,phi)
pow_ex_nr <- POW_EX_NR(m,n,kx_y,CritC,p1,p2)
POWEX <- pow_ex_nr[1,1] 
if (POWEX >= POWO) break                      }
  
  m <- m-incr5_m
  n <- n-incr5_n
  
  repeat
  { m <- m+incr_m
  n <- n+incr_n
  CritC <- CritC_condS(m,n,alpha)
  phi <- CritfctXY(m,n,CritC)
  kx_y <- KX_Y(m,n,phi)
  pow_ex_nr <- POW_EX_NR(m,n,kx_y,CritC,p1,p2)
  POWEX <- pow_ex_nr[1,1]
  if (POWEX >= POWO) break                }   }  else
  { repeat 
  { m <- m-incr5_m
  n <- n-incr5_n
  CritC <- CritC_condS(m,n,alpha)
  phi <- CritfctXY(m,n,CritC)
  kx_y <- KX_Y(m,n,phi)
  pow_ex_nr <- POW_EX_NR(m,n,kx_y,CritC,p1,p2)
  POWEX <- pow_ex_nr[1,1]
  if (POWEX < POWO) break                      }
    
    repeat
    { m <- m+incr_m
    n <- n+incr_n
    CritC <- CritC_condS(m,n,alpha)
    phi <- CritfctXY(m,n,CritC)
    kx_y <- KX_Y(m,n,phi)
    pow_ex_nr <- POW_EX_NR(m,n,kx_y,CritC,p1,p2)
    POWEX <- pow_ex_nr[1,1]
    if (POWEX >= POWO) break                       }   }

Mex <- m
Nex <- n


POWNR <- pow_ex_nr[1,2]
mstart <- Mex
nstart <- Nex

if (POWNR < POWO)
{ repeat
{ m <- m+incr5_m
n <- n+incr5_n
CritC <- CritC_condS(m,n,alpha)
phi <- CritfctXY(m,n,CritC)
kx_y <- KX_Y(m,n,phi)
pow_ex_nr <- POW_EX_NR(m,n,kx_y,CritC,p1,p2)
POWNR <- pow_ex_nr[1,2]
if (POWNR >= POWO) break                      }
  m <- m-incr5_m
  n <- n-incr5_n
  repeat
  { m <- m+incr_m
  n <- n+incr_n
  CritC <- CritC_condS(m,n,alpha)
  phi <- CritfctXY(m,n,CritC)
  kx_y <- KX_Y(m,n,phi)
  pow_ex_nr <- POW_EX_NR(m,n,kx_y,CritC,p1,p2)
  POWNR <- pow_ex_nr[1,2]
  if (POWNR >= POWO) break                }   }  else
  { repeat
  { m <- m-incr5_m
  n <- n-incr5_n
  CritC <- CritC_condS(m,n,alpha)
  phi <- CritfctXY(m,n,CritC)
  kx_y <- KX_Y(m,n,phi)
  pow_ex_nr <- POW_EX_NR(m,n,kx_y,CritC,p1,p2)
  POWNR <- pow_ex_nr[1,2]
  if (POWNR < POWO) break                       }
    repeat
    { m <- m+incr_m
    n <- n+incr_n
    CritC <- CritC_condS(m,n,alpha)
    phi <- CritfctXY(m,n,CritC)
    kx_y <- KX_Y(m,n,phi)
    pow_ex_nr <- POW_EX_NR(m,n,kx_y,CritC,p1,p2)
    POWNR <- pow_ex_nr[1,2] 
    if (POWNR >= POWO) break                     }   }

Mnr <- m
Nnr <- n


m0 <- (Mex+Mnr)/2
n0 <- (Nex+Nnr)/2
midN_m <- round((m0+n0)/mton_ab + 10**-6)*mton_a
midN_n <- round((m0+n0)/mton_ab + 10**-6)*mton_b


erg <- c(s_mstart,s_nstart,Mex,Nex,POWEX,Mnr,Nnr,POWNR,midN_m,midN_n)

}
