\name{siland}
\alias{siland}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of spatial influence of landscape
}
\description{
siland is used to fit spatial influence of landscape.  
%It can be used to carry out gaussian, exponential or binomial models. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
siland(loc.model, land=NULL, data, initSIF = NULL, sif = "exponential", 
family = "gaussian", test = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{loc.model}{
a symbolic description (see lm() or glm()) of the response variable concerning local variables. Random effects are also allowed according to the syntax in package lme4 (see lmer() function in package lme4).
}
  \item{land}{
a list containing the location of the landscape variables. Each component of the list is a dataframe with two columns "X" and "Y" indicating the locations of pixels for a given landscape variable.
}
  \item{data}{
a dataframe containing the response variable and the local variables.
}
  \item{initSIF}{
a vector indicating the starting values for the estimation of the mean distance of the spatial influence functions. The length of initSIF has to be equal to the length of list land. By default, initialisation is equal to 100 for each landscape variable.
}
  \item{sif}{
the family of the spatial influence function. sif can be "exponential" or "gaussian". 
}
  \item{family}{
the distribution of response variable. family can be "gaussian", "poisson" or "binomial" and the associated link function are identity, log and logit respectively.
}
  \item{test}{
logical value. If TRUE, a ratio likelihood test is performed for each explanatory variable (local and spatial). By default, the value is FALSE since performing tests for all the variables can be heavy computing.
}
%%  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
%%}
}

\value{
siland returns an object of type list.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{coefficients}{vector of estimated coefficients}
\item{local}{an object of class formula that indicates the local model used}
\item{landcontri}{a dataframe of estimated contributions of each spatial variable (in column) to each observation (in row). The number of columns is equal to the length of list land}  
\item{loglik}{log-likelihood for the estimated parameters}     
\item{loglik0}{log-likelihood for the local model}      
\item{fitted}{fitted values}
\item{sif}{the family of the spatial influence function}          
\item{resoptim}{an object of class optim or optimize giving informations about the optimization procedure see optim() or optimize() for further details.}     
\item{AIC}{akaike information criterion}         
\item{AIC0}{akaike information criterion for local model (no landscape variable)}         
\item{nparam}{number of parameters }       
\item{pval0}{p.value of the test of global effect of spatial variables. Obtained from the likelihood ratio test between the complete model and the local model.}   
\item{pval}{if test is TRUE, vector of p.values of the test of effect of each landscape variables}    
\item{family}{family distribution for the model}
\item{sd.error}{standard error for gaussian family, NA in other case}
\item{model.Type}{type of local model: GLM for generalised model, LMM for linear mixed model or GLMM for generalised linear mixed model}
\item{rand.StdDev}{standard deviation of random effects for LMM or GLMM}
\item{nparam}{number of model parameters}
\item{err}{estimated residuals}

}

\author{
Carpentier, F. and Martin, O.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\references{
Chandler R. and Hepinstall-Cymerman J. (2016) Estimating the spatial scales of landscape effects on abundance. Landscape ecology, 31: 1383-1394.
   }

\examples{
\donttest{
data(dataSiland)
data(landSiland)
resE=siland(y~locvar,land=landSiland,data=dataSiland,sif="exponential",family=gaussian)
resE
resE$AIC
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }%
