\name{Bsiland}
\alias{Bsiland}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of landscape influence with buffers
}
\description{
Bsiland is used to find the buffer size to estimate lanbdscape influence.
}
\usage{
Bsiland(formula, land, data, family = "gaussian", initb = 50, border = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a symbolic description (see lm() or glm()) of the response variable concerning local variables. Random effects are also allowed according to the syntax in package lme4 (see lmer() function in package lme4).
}
  \item{land}{
an object of class sf that gives the landscape variables for the model.
}
  \item{data}{
a dataframe containing the response variable and the local variables.
}
  \item{family}{
the distribution of response variable. family can be "gaussian", "poisson" or "binomial" and the associated link function are identity, log and logit respectively.
}
  \item{initb}{
a vector indicating the starting values to estimate the buffer sizes. By default, buffer sizes are initialized to 100 for each variable.
}
  \item{border}{
a logical indicating whereas the buffer are computed from the observation locations (border=F) or from the border of the polygon where observations are located (border=T)}
}


\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Fsiland returns an object of type list.

\item{coefficients}{vector of estimated coefficients}
\item{parambuffer}{vector of estimated buffer distances}
\item{formula}{an object of class formula that indicates the local model used}
\item{buffer}{a dataframe with percentage for each landscape variable inside buffers}
\item{loglik}{log-likelihood for the estimated parameters}     
\item{loglik0}{log-likelihood for the local model} 
\item{fitted}{fitted values}
\item{resoptim}{an object of class optim or optimize giving informations about the optimization procedure see optim() or optimize() for further details.}     
\item{result}{an object of type lm/glm/lmer that corresponds to the estimated model conditionnaly to the best buffer sizes for the landscape variables.}
\item{AIC}{akaike information criterion}         
\item{AIC0}{akaike information criterion for local model (no landscape variable)}         
\item{nparam}{number of parameters }       
\item{pval0}{p.value of the test of global effect of spatial variables. Obtained from the likelihood ratio test between the complete model and the local model.}   
\item{family}{family distribution for the model}
\item{sd.error}{standard error for gaussian family, NA in other case}
\item{model.Type}{type of local model: GLM for generalised model, LMM for linear mixed model or GLMM for generalised linear mixed model}
\item{rand.StdDev}{standard deviation of random effects for LMM or GLMM}
\item{err}{estimated residuals}
\item{newdata}{a dataframe with the local variables and the percentages for the different landscape variables obtained with the best buffer sizes.}
\item{border}{a logical indicating the value used for estimation}

}


\author{
O. Martin and F. Carpentier
}

\examples{
\donttest{
data(dataSiland)
data(landSiland)
resB=Bsiland(obs~x1+L1+L2,land=landSiland,data=dataSiland,initb = c(50))
resB
summary(resB)
Bsiland.lik(resB,land=landSiland,data=dataSiland,varnames=c("L1","L2"),seqd=seq(50,1000,length=20))


}

}

%% ~Make other sections like Warning with \section{Warning }{....} ~



