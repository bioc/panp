"pa.calls" <-
function(object=NULL, looseCutoff=0.02, tightCutoff=0.01, verbose=FALSE) {

# Print intro
if (is.null(object)) {
    cat("\nUSAGE:
\tpa.calls(object, looseCutoff=0.02, tightCutoff=0.01, verbose = FALSE)\n
INPUTS: 
\tobject - exprSet; e.g., returned from expression-generating function, 
\t   such as expresso(), rma(), mas5(), etc. 
\tlooseCutoff - the larger P-value cutoff
\ttightCutoff - the smaller, more strict P-value cutoff
\tverbose - TRUE or FALSE\n
OUTPUTS: 
\tReturns a new exprSet, with exprs = Pcalls, se.exprs = Pvals:
\tPvals - a matrix of P-values of same dimensions as exprs(exprSet). Each
\t   datapoint is the P-value for the probeset at the same x,y coordinates. 
\tPcalls - a matrix of Presence (P), Marginal (M), Absent (A) indicators\n\n")
return()
}

# If we have inputs, let's continue:

# First, make sure that "affy" library is installed
if (require(affy)==FALSE){
  stop("\npa.calls() requires BioConductor 'affy' package.
  Currently, it is not installed. Please install and load, then
  rerun pa.calls()\n\n")
}
if(verbose){
   cat("\nInvoking function 'pa.calls'\n")
   cat("tightCutoff is ",tightCutoff,"\nlooseCutoff is ",looseCutoff,"\n")
}
# check inputs for correct type and range
if (class(object)[1]!="exprSet") {
stop("\nAborting: object must be an exprSet\n\n")
}
# exprSet type MUST be only hgu133a or hgu133plus2; no other chipset
# is currently supported
chip = annotation(object)
if (chip == "hgu133b"){
stop("\nHG-U133B is not currently supported. Supported chip types are\n
HG-U133A and HG-U133 Plus 2.0\n\n")
}
if ((chip != "hgu133a")&(chip != "hgu133atag")&(chip != "hgu133plus2")){
 stop("\nAborting: chip type must be either HG-U133A or HG-U133 Plus 2.0 \n\n")
}
# Check cutoff ranges and relationship
if ((tightCutoff > 1.0)|(tightCutoff < 0.0)|
(looseCutoff > 1.0)|(looseCutoff < 0.0)) {
 stop("\nAborting: cutoffs must be between 0.0 and 1.0\n\n")
}
if (tightCutoff > looseCutoff){ 
 stop("\nAborting: tightCutoff must be lower than looseCutoff\n\n")
}
if(verbose){
  cat("Outputs will be exprSet containing two matrices of same dimensions as input full dataset:
  1. se.exprs = Pvals - a matrix of P-values 
  2. exprs = Pcalls - a matrix of P/A/M indicators: 
  'P': P-values <= tightCutoff ",tightCutoff,"
  'A': P-values > looseCutoff ",looseCutoff,"
  'M': P-values between ", tightCutoff, " and ", looseCutoff,"\n")
}
# Extract exprs and negative exprs
# First, read in the correct chip-specific NSMP list
if ((chip == "hgu133a")|(chip == "hgu133atag")) {
data(NSMPnames.hgu133a)
NSMPnames <- NSMPnames.hgu133a
rm(NSMPnames.hgu133a, envir = globalenv())
}else{
data(NSMPnames.hgu133plus2)
NSMPnames <- NSMPnames.hgu133plus2
rm(NSMPnames.hgu133plus2, envir = globalenv())
}
# Now the object "NSMPnames" has the chip-specific list of NSMPs
# Get all exprs, then get just the NSMP subset
AllExprs <- exprs(object)
NegExprs <- AllExprs[NSMPnames,]
AllExprs<- as.matrix(AllExprs)
NegExprs<-as.matrix(NegExprs)

# Calculations…

# Create cutoff function
cutoff_fcn <- function(x){
 if (x<=tightCutoff) return("P") 
 else if (x > looseCutoff) return("A")
else return("M")
}

# set up loop to do each chip in exprSet separately
len <- length(colnames(AllExprs))
cat("\nProcessing",len, "chips: ")
flush.console()

for (chip in 1:len) {
xNeg <- NegExprs[,chip]
xNeg <- sort(xNeg, decreasing = TRUE)# sort the NSMP exprs, decreasing
yNeg <- seq(0,1,1/(length(xNeg)-1))# generate the y axis for the survivor fcn.
interp <- approxfun(xNeg, yNeg, yleft=1, yright=0)# create interpolation fcn.
PV <- sapply(AllExprs[,chip], interp)# calculate column of P-values
PC <- sapply(PV, cutoff_fcn)# create P/A/M column
if (chip==1) {
Pvals <- PV
Pcalls <- PC
 }else{
Pvals <- cbind(Pvals, PV)
Pcalls <- cbind(Pcalls, PC)
}
cat("#")
flush.console()
} # END LOOP
Pvals <- as.matrix(Pvals)
Pcalls <- as.matrix(Pcalls)
colnames(Pvals) <- colnames(AllExprs)
colnames(Pcalls) <- colnames(AllExprs)
eset <- new("exprSet",
              exprs = Pcalls,
              se.exprs = Pvals,
              phenoData = phenoData(object),
              description = description(object),
              annotation = annotation(object),
              notes = c(notes(object)))
cat("\nProcessing complete.\n\n")
flush.console()

# Print out intensity values at the two cutoffs:
myX <- NegExprs
maxLen <- 0
for (i in 1:len) {
stringLen <- nchar(colnames(AllExprs)[i])  # length of column name in chars
if (stringLen > maxLen) maxLen <- stringLen
}
maxLen <- maxLen-6
if(maxLen<2) {maxLen=2}
for (i in 1:len) {
myY <- seq(0,1,1/(length(myX[,i])-1))
myX[,i] <- sort(myX[,i], decreasing = TRUE)
revInterp <- approxfun(myY, myX[,i], yleft=1, yright=0)
revTight <- revInterp(tightCutoff)
revLoose <- revInterp(looseCutoff)
if(i==1) {
cat("\nIntensities at cutoff P-values of ", looseCutoff," and ", tightCutoff, ":\n")
cat("Array:")
for(j in 1:maxLen) { cat(" ") }
cat("value at",looseCutoff,"   value at",tightCutoff,"\n")
} 
#...and for all i:
cat(colnames(AllExprs)[i], "\t", format.pval(revLoose,digits=3), "\t\t", format.pval(revTight,digits=3), "\n")
}
cat("\n")
cat("[NOTE: 'Collapsing to unique x values...' warning messages are benign.]\n\n")
return(eset)  # Pcalls are in eset$exprs slot ; Pvals are in eset$se.exprs slot
}

