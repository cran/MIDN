\name{fisher_boschloo_midN}
\alias{fisher_boschloo_midN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nearly exact sample size calculation for the Fisher-Boschloo test for differences between 
independent binomial proportions
}
\description{
The function computes the exact sample sizes required in the randomized UMPU test and its conservative nonrandomized 
version for attaining prespecified power. In a final step, the mean of both quantities is output as an nearly exact value required 
in the Fisher-Boschloo test, a powerful nonrandomized version of the exact Fisher-type test.
}
\usage{
fisher_boschloo_midN(alpha, SW, p1, p2, POWO, mton_a, mton_b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{
target significance level
}
  \item{SW}{
step width for increasing p2 in the search for the size of     
a given critical region in the sample space of (X,Y)
}
  \item{p1}{
true value of the responder rate for Population 1
}
  \item{p2}{
true value of the responder rate for Population 2
}
  \item{POWO}{
power to be obtained against the alternative (p1,p2)
}
  \item{mton_a}{
desired ratio of sample sizes: numerator
}
  \item{mton_b}{
desired ratio of sample sizes: denominator
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{mstart}{initial value of 1st sample size}
  \item{nstart}{initial value of 2nd sample size}
  \item{Mex}{size of Sample 1 for randomized UMPU test}
  \item{Nex}{size of Sample 2 for randomized UMPU test}
  \item{POWEX}{power of randomized UMPU test attained with m=Mex,n=Nex}
  \item{Mnr}{size of Sample 1 for conservative nonrandomized Fisher-type test}
  \item{Nnr}{size of Sample 2 for conservative nonrandomized Fisher-type test}
  \item{POWNR}{power of conservative nonrandomized Fisher-type test attained with 
        m=Mnr,n=Nnr}
  \item{midN_m}{nearly exact size of Sample 1 for Boschloo-Fisher test}
  \item{midN_n}{nearly exact size of Sample 1 for Boschloo-Fisher test}
}
\references{
Wellek S: Nearly exact sample size calculation for powerful nonrandomized tests for differences between 
binomial proportions. Statistica Neerlandica 69 (2015), 358-373.
}
\author{
Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

result1 <- fisher_boschloo_midN(0.025,0.0001,0.95,0.8,0.8,2,1)
POWEX <- result1[5]
result1  # shows values of vector result1
POWEX    # shows value of POWEX

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ binomial two-sample problem }
\keyword{ Boschloo's approach }
\keyword{ exact Fisher-type test }

