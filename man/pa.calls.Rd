\name{pa.calls}
\alias{pa.calls}
\title{Presence-Absence Calls from Negative Strand Matching Probesets}
\description{
Function to make gene presence/absence calls based on distance from empirical distribution of chip-specific negative strand matching probesets (NSMP).
}
\usage{
pa.calls(object, looseCutoff = 0.02, tightCutoff = 0.01, verbose = FALSE)
}
\arguments{
  \item{object}{an exprSet object (result of running expression-generating function, like expresso(), rma(), mas5(), etc.) Currently, this must be of chip type HGU133A or HGU133 Plus 2.0}
  \item{looseCutoff}{the larger P-value cutoff (see details)}
  \item{tightCutoff}{the smaller, more strict P-value cutoff}
  \item{verbose}{logical. If 'TRUE' detailed progress messages are reported.}
}
\details{The function calculates a matrix of P-values for the expression values in the input exprSet. P-values are calculated based on the empirical survivor function (1-CDF) of the set of negative probesets identified by Affymetrix as negative strand matching probesets (NSMP) with no cross hybridization. These probesets are therefore assumed to show nothing but background/machine noise plus some occasional non-specific binding. The P-value returned for any probeset expression value in exprSet is the value of the NSMP survivor function for that expression level. 

Presence/Absence calls are derived by applying the two cutoff values to
the matrix of P-values for all genes in the exprSet, as follows:
\item Present ('P'):	P-values <= tightCutoff
\item Absent ('A'):	P-values > looseCutoff
\item Marginal ('M'):	P-values between tightCutoff and looseCutoff	
}
\value{
  \item{exprSet}{a new exprSet: exprs slot holds Pcalls; se.exprs slot holds Pvals, as follows:}
  \item{Pvals }{a matrix of P-values. Each data point is the P-value for the expr 
at the same x,y coordinates. }
  \item{Pcalls }{a matrix of Presence (P), Marginal (M), Absent (A) indicators}
}
\note{NSMP sets have been established for the HGU133A and HGU133-Plus-2.0 chipsets to date. Hence only these two are currently supported by PANP.}
\references{Warren, P., Bienkowska, J., Martini, P., Jackson, J., and Taylor, D., PANP - a New Method of Gene Detection on Oligonucleotide Expression Arrays (2005), in preparation}
\author{Peter Warren}

\examples{
## Load example exprSet 
data(test133a)

## Generate Pvals and Pcalls matrices from exprSet, using default cutoffs
newExprSet <- pa.calls(test133a)

## to access the Pcalls and Pvals:
myPcalls <- exprs(newExprSet)
myPvals <- se.exprs(newExprSet)
}
\keyword{manip}

