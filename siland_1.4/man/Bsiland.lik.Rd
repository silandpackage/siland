\name{Bsiland.lik}
\alias{Bsiland.lik}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute -(Log-likelihood) for models with different buffer sizes for landscape variables.
}
\description{
This function allows to investigate some optimization problems for function Bsiland by plotting the -(log-likelihhod) 
for models with different buffer sizes for landscape variables.
}
\usage{
Bsiland.lik(res, land,data, varnames = NULL, seqd = seq(2, 2000, length = 10))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{
an object from Bsiland estimation
}

  \item{land}{
the landscape (a sf object) used to obtain the result res.
}

  \item{data}{
the dataframe used to obtain the result res.
}
  \item{varnames}{
landscape variables for which -(log-likelihhod) is studied. If varnames is null (by default), -(log-likelihhod) is investigated for 
all landscape variables in estimted model for object res.
}
  \item{seqd}{
a vector of positives values that represent different buffer sizes for the landscape variables.
}
}
\details{
Estimation based on log-likelihood aims to find parameters that maximises log-likelihood, that is to say 
that minimises -(log-likelihood). From graphics obtained with this function, orange curve gives 
the -(log-likelihood) obtained for the estimated model. Others curves gives -(log-likelihood) by varying 
buffer sizes for landscape variables. If minimisation for -(log-likelihood) obtained for object res occured correctly, no curve should go below the orange one. Otherwise, the obtained minimisation corresponds to a local minimum, and a new estimation should be proceed by changing starting point (see argument init in function Bsiland).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

A matrix with -(log-likelihood) values for the different landscape variables (argument varnames) and for the different buffer sizes (argument seqd).
}

\author{
O. Martin and F. Carpentier
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
\donttest{
data(dataSiland)
data(landSiland)
resB=Bsiland(obs~x1+L1+L2,land=landSiland,data=dataSiland,initb = c(50))
Bsiland.lik(resB,dataSiland,land=landSiland,varnames=c("L1","L2"),seqd=seq(5,500,length=20))
}
}



