\name{MIDN-package}
\alias{MIDN-package}
\alias{MIDN}
\docType{package}
\title{
Nearly exact sample size calculation for exact powerful nonrandomized tests
       for differences between binomial proportions
}
\description{
Implementation of the mid-n algorithms presented in 
             Wellek S (2015) Statistica Neerlandica 69, 358-373 for exact 
             sample size calculation for superiority trials with binary outcome.
}
\author{
Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de> \cr \cr
Maintainer: \cr
            Peter Ziegler <peter.ziegler@zi-mannheim.de> \cr
}
\references{
Wellek S: Nearly exact sample size calculation for powerful nonrandomized tests for differences between 
binomial proportions. Statistica Neerlandica 69 (2015), 358-373.
}
\examples{
result1 <- fisher_boschloo_midN(0.025,0.0001,0.95,0.8,0.8,2,1)
POWEX <- result1[5]
result1  # shows values of vector result1
POWEX    # shows value of POWEX

result2 <- McNem_Score_midn(0.025,0.0001,0.585,0.315,0.9)
POWEX <- result2[3]
result2  # shows values of vector result2
POWEX    # shows value of POWEX
}
