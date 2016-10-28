\name{McNem_Score_midn}
\alias{McNem_Score_midn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nearly exact sample size calculation for the level-corrected score test for differences 
between binomial proportions estimated from paired data
}
\description{
Again, the function computes the exact sample sizes required in the 
randomized UMPU test and its conservative nonrandomized counterepart 
for attaining prespecified power. However, in contrast to the parallel 
group setting, the midpoint of the interval between these two numbers 
shall now used as an nearly exact value of the number of pairs to be observed 
in the asymptotic test based on the score-statistic corrected for possible 
exceedances of the nominal significance level.
}
\usage{
McNem_Score_midn(alpha, SW, ppl, pmi, POWO)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{
target significance level, 1-sided
}
  \item{SW}{
width of search grid for determining the size of a given critical
region in the sample space of N+ [= number of pairs with (Xi,Yi) = (1,0)]
and N0 [= number of tied pairs]
}
  \item{ppl}{
true value of Pr[(X,Y) = (1,0)]
}
  \item{pmi}{
true value of Pr[(X,Y) = (0,1)]
}
  \item{POWO}{
power to be attained in the level-corrected score test against the 
alternative (ppl,pmi)
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{ nstart }{initial value for the iteration algorithm}
  \item{ Nex }{sample size required in the exact randomized McNemar test}
  \item{ POWEX }{power of the exact randomized McNemar test performed with Nex pairs}
  \item{ Nnr }{sample size required in the conservative nonrandomized McNemar test}
  \item{ POWNR }{power of the nonrandomized McNemar test performed with Nnr pairs}
  \item{ mid_n }{midpoint of the interval [Nex,Nnr], rounded to the next integer}
}
\references{
Wellek S: Nearly exact sample size calculation for powerful nonrandomized tests for differences 
between binomial proportions. Statistica Neerlandica 69 (2015), 358-373.
}
\author{
Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>
}

\examples{

result2 <- McNem_Score_midn(0.025,0.0001,0.585,0.315,0.9)
POWEX <- result2[3]
result2  # shows values of vector result2
POWEX    # shows value of POWEX

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ exact nonconditional test}
\keyword{ McNemar setting}
\keyword{ score statistic }
