\name{Fsiland}
\alias{Fsiland}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of landscape influence with spatial influnce functions
}
\description{
Fsiland is used to fit spatial influence of landscape.  
%It can be used to carry out gaussian, exponential or binomial models. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Fsiland(formula,land,data,family ="gaussian",sif="exponential", init = 100,
 border=F,wd=50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a symbolic description (see lm() or glm()) of the response variable according local and landscape variables. Random effects are also allowed according to the syntax in package lme4 (see lmer() function in package lme4).
}
  \item{land}{
an object of class sf that gives the landscape variables for the model
}
  \item{data}{
a dataframe containing the response variable and the local variables.
}

\item{family}{
the distribution of response variable. family can be "gaussian", "poisson" or "binomial" and the associated link function are identity, log and logit respectively.
}
\item{sif}{
the family of the spatial influence function. sif can be "exponential", "gaussian" or "uniform".
}
  \item{init}{
a vector indicating the starting values for the estimation of the mean distance of the spatial influence functions. If one value is given, parameter for each spatial influnce function are initialised with this value. By default, initialisation is equal to 100 for each landscape variable.
}
  
  \item{border}{
a logical indicating if pixels into the field (polygon) where obervations are located have to be take into account to evaluate the influences. If border=T (by default), all pixels are take into account. If border=F, only pixels outside the fields (polygon) are condidered to evaluate influence for landscape variables.
}
  
\item{wd}{
a number indiacting the size of pixels for raster discretisazion.
}
  
%%  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
%%}
}

\value{
Fsiland returns an object of type list.

\item{coefficients}{vector of estimated coefficients}
\item{paramSIF}{vector of estimated parameters for spatial influence functions}
\item{formula}{an object of class formula that indicates the model}
\item{landcontri}{a dataframe of estimated contributions of each spatial variable (in column) to each observation (in row).}  
\item{loglik}{log-likelihood for the model}     
\item{loglik0}{log-likelihood for the local model (only local variables)} 
\item{result}{an object of type lm/glm/lmer that corresponds to the estimated model conditionnaly to the parameters for the spatial influence functions. } 
\item{fitted}{fitted values}
\item{sif}{the family of the spatial influence function}          
\item{resoptim}{an object of class optim or optimize giving informations about the optimization procedure see optim() or optimize() for further details.}     
\item{AIC}{akaike information criterion}         
\item{AIC0}{akaike information criterion for local model (no landscape variable)}         
\item{nparam}{number of parameters }       
\item{pval0}{p.value of the test of global effect of spatial variables. Obtained from the likelihood ratio test between the complete model and the local model.}   
\item{family}{family distribution for the model}
\item{sd.error}{standard error for gaussian family, NA in other case}
\item{model.Type}{type of local model: GLM for generalised model, LMM for linear mixed model or GLMM for generalised linear mixed model}
\item{rand.StdDev}{standard deviation of random effects for LMM or GLMM}
\item{err}{estimated residuals}
\item{border}{a logical indicating the value used for estimation}
\item{wd}{a number indicating the value used for raster discretisazion}
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
resF=Fsiland(obs~x1+L1+L2,land=landSiland,data=dataSiland,sif="exponential")
resF
resF$AIC
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }%
