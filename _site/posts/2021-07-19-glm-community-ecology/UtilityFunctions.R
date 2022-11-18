#########################################################
## Utility function to compute the different correlations
#########################################################
library(ade4)
TraitEnvCor <- function(L,E,T, Chessel = TRUE){
  
  E<-as.matrix(E)
  T<-as.matrix(T)
  L<-as.matrix(L)
  #  centering_mat <- function(X,w){ X - rep(1,length(w))%*%t(w)%*% X }
  standardize_w <- function(X,w){
    ones <- rep(1,length(w))
    Xc <- X - ones %*% t(w)%*% X
    Xc / ones%*%sqrt(t(ones)%*%(Xc*Xc*w)) 
  } 
  
  # check_L()
  rows<-seq_len(nrow(L))
  cols<-seq_len(ncol(L))
  rni <-which(rowSums(L)==0)
  repeat {
    if (length(rni)) {L <- L[-rni,,drop = FALSE]; rows <-rows[-rni]}
    ksi <- which(colSums(L)==0)
    if (length(ksi)) {L <- L[,-ksi, drop = FALSE]; cols <- cols[-ksi]}
    rni <-which(rowSums(L)==0)
    if ( length(rni)==0 & length(ksi)==0){break}
  }
  E <-E[rows,,drop = FALSE]
  T <-T[cols,,drop = FALSE]
  # end check_L()
  
  L<-L/sum(L)
  # dimensions
  #S <- ncol(L) # number of species
  #n <- nrow(L) # number of communities
  p <- ncol(E) # number of environmental predictors
  q <- ncol(T) # number of traits
  
  # setting up matrices
  Wn <- rowSums(L)
  Ws <- colSums(L)
  # cor matrices are trait by environment
  CWM <- L%*%T/Wn  # weighted means wrt to T
  CWM.cor <- cor(CWM,E)
  
  SNC <- t(L)%*%E/Ws  # weighted means wrt to E
  SNC.cor <- cor(T,SNC)
  
  CWMstd_w  <- standardize_w(CWM,Wn)
  Estd_w <- standardize_w(E,Wn)
  wCWM.cor <- t(t(Estd_w)%*%(CWMstd_w*Wn))
  
  SNCstd_w <- standardize_w(SNC,Ws)
  Tstd_w <-  standardize_w(T,Ws)
  wSNC.cor <- t(Tstd_w)%*%(SNCstd_w*Ws)
  
  # Fourth corner calculated as W_n weighted covariance between 
  # CWM and standardized T (trait)
  
  CWM_std_tw <- L%*%Tstd_w/Wn #CWM wrt to standardized T (trait)
  Fourthcorner <- t(CWM_std_tw)%*%(Estd_w*Wn)
  if (Chessel){
    singular_val1 <- sqrt(ade4::dudi.coa(L, scannf = FALSE)$eig[1])
    Chessel.4thcor<-Fourthcorner/ singular_val1
  }else { Chessel.4thcor<-NA;singular_val1 <-1}
  
  
  # variation components
  # Among communities 
  Among.Variation <- sum(diag(t(CWM_std_tw)%*%(CWM_std_tw* Wn)))
  # Within communities 
  Within.Variation <- 1 - Among.Variation
  
  # result specialized to one trait and one environment variables; use array(0, dim(6,k,p)) in the general case
  # array.result<-matrix(c(CWM.cor,wCWM.cor,SNC.cor,wSNC.cor,Fourthcorner,Chessel.4thcor,Mean.Variation),ncol=1)
  array.result<-array(0, dim=c(8,q,p))
  rownames(array.result)<- c("CWM.cor","wCWM.cor","SNC.cor","wSNC.cor","Fourthcorner","Chessel.4thcor","Among Wn-variance (%)", "Within Wn-variance (%)")
  array.result[1,,]<-CWM.cor
  array.result[2,,]<-wCWM.cor
  array.result[3,,]<-SNC.cor
  array.result[4,,]<-wSNC.cor
  array.result[5,,]<-Fourthcorner
  array.result[6,,]<-Chessel.4thcor
  array.result[7,,]<-Among.Variation * 100
  array.result[8,,]<-Within.Variation * 100
  return(array.result[,,])
}

###################################################################
## Utility function for the row, column, and row-column permutation schemes
## for a single trait and a single environmental variable
## and the five test statistics/approaches of the paper
## CWM.cor, wCWM.cor, SNC.cor wSNC.cor and Fourthcorner
###################################################################

CorPermutationTest <- function(L, E, T, nrepet = 999){
  E<-as.matrix(E)
  T<-as.matrix(T)
  L<-as.matrix(L)
  obs <- TraitEnvCor(L,E,T)[1:5]
  sim.row <- matrix(0, nrow = nrepet, ncol = ncol(E) * 5)
  sim.col <- matrix(0, nrow = nrepet, ncol = ncol(E) * 5)
  for(i in 1:nrepet){
    per.row <- sample(nrow(L))
    per.col <- sample(ncol(L))
    sim.row[i, ] <- c(as.matrix(data.frame(TraitEnvCor(L,E[per.row,,drop= FALSE],T))))[1:5]
    sim.col[i, ] <- c(as.matrix(data.frame(TraitEnvCor(L,E,T[per.col,,drop= FALSE]))))[1:5]
  }
  pval.row <- (rowSums(apply(sim.row^2, 1, function(i) i >= obs^2)) + 1)  / (nrepet + 1)
  pval.col <- (rowSums(apply(sim.col^2, 1, function(i) i >= obs^2)) + 1)  / (nrepet + 1)
  
  result <- cbind(cor = obs, prow = pval.row, pcol = pval.col, pmax = apply(cbind(pval.row, pval.col), 1, max))
  return(result)
}

