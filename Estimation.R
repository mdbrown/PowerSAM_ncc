








getEstandSE_SP <- function(
                        data, 
                       cutpoint,
                       predict.time, 
                        nmatch)
{  
  
#  browser()
  
 # N = nrow(data)
  #data$vi = 1; data$wi = 1

  #decompose data into its three parts
  mydata.ncc = data[[1]];
  Vij.IND = data[[2]]; 
  Iik0.mat = data[[3]]
  B0 = 200
  junk=NCC.Cox.explicit.PTB.s.y(data = mydata.ncc,
                                V.IND = Vij.IND,
                                Iik0 = Iik0.mat,
                                wgtk.ptb = NULL,
                                B0 = B0,
                                Zmatch.ind = NULL, 
                                t0 = predict.time, 
                                m.match = nmatch);
  
  Ptb.result = matrix(0,B0,5)
  
  for (i in 1:B0) {
    Ptb.result[i,] <- unlist(Ptb.ROC.FUN(junk$ck,
                                         junk$CondSck.ptb[,i],
                                         junk$Fck.ptb[,i],
                                         uu0 = cutpoint,
                                         type = "cutoff", 
                                         junk$yk))
  }
  
  
  est <- c(junk$beta.est, colMeans(Ptb.result, na.rm = TRUE))
  se  <- c(junk$beta.sd,apply(Ptb.result,2,sd))

  est <- as.data.frame(t(est)); names(est) =c("coef", "AUC", "FPR","TPR","NPV","PPV")
  se <- as.data.frame(t(se)); names(se) =c("coef", "AUC", "FPR","TPR", "NPV","PPV")
  
  list(estimates = unlist(est),se =se) 
     
}



getEstandSE_NP <- function(
  data, 
  cutpoint,
  predict.time, 
  nmatch)
{  
  
    browser()
  
  # N = nrow(data)
  #data$vi = 1; data$wi = 1
  
  #decompose data into its three parts
  mydata.ncc = data[[1]];
  Vij.IND = data[[2]]; 
  Iik0.mat = data[[3]]
  B0 = 200
  junk=NCC.NP.kn.PTB.s.y(data = mydata.ncc,
                                V.IND = Vij.IND,
                                Iik0 = Iik0.mat,
                                wgtk.ptb = NULL,
                                B0 = B0,
                                Zmatch.ind = NULL, 
                                t0 = predict.time,
                                m.match=nmatch);
  
  Ptb.result = matrix(0,B0,5)
  
  for (i in 1:B0) {
    Ptb.result[i,] <- unlist(Ptb.ROC.FUN(junk$ck,
                                         junk$CondSck.ptb[,i],
                                         junk$Fck.ptb[,i],
                                         uu0 = cutpoint,
                                         type = "cutoff", 
                                         junk$yk))
  }
  
  
  est <- c(junk$beta.est, colMeans(Ptb.result, na.rm = TRUE))
  se  <- c(junk$beta.sd,apply(Ptb.result,2,sd))
  
  est <- as.data.frame(t(est)); names(est) =c("coef", "AUC", "FPR","TPR","NPV","PPV")
  se <- as.data.frame(t(se)); names(se) =c("coef", "AUC", "FPR","TPR", "NPV","PPV")
  
  list(estimates = unlist(est),se =se) 
  
}




VTM <- function(vc, dm)
{
    matrix(vc, ncol = length(vc), nrow = dm, byrow = T)
}

