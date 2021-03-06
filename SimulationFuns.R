
PowerSim <- function(N, parameter, S.0, t.0,
                             cutoff=NULL, 
                             predict.time,
                             parval.H0, 
                             parval.Ha, 
                             ESTmethod, 
                             f_x = dnorm, 
                             F_xInv = qnorm, 
                             mm = 5,  
                             alpha = 0.05, 
                             cens.lam = .2, 
                             time.max = NULL,
                              nmatch = 3){
   a = -log(S.0)/t.0

  mybetas <- get.Betas(parameter, a, predict.time, cutoff, parval.H0, parval.Ha)
  beta.H0 <- mybetas[1]
  beta.Ha <- mybetas[2]
   if(ESTmethod == "SP"){
     mynames <- c("beta", "AUC", "TPR", "FPR" , "PPV", "NPV")
   }else{
     mynames <- c( "AUC", "TPR", "FPR" , "PPV", "NPV")
   }
   
   EST <- matrix(nrow = mm, ncol = length(mynames))
 
   SE_ncc <- matrix(nrow = mm, ncol = length(mynames))
  
  N_ncc <- numeric(mm)
   for(i in 1:mm){

    tmpDat <- SIM.data.singleMarker(nn = N, 
                                    beta = beta.Ha, 
                                    lam0 = a, 
                                    cens.lam = cens.lam, 
                                    time.max = time.max, 
                                    m.match = nmatch) 
    tmpDat$data$id <- 1:N
    
    estimates_ncc <- survMTP.ncc(time = xi,
                                 event = di, 
                                 marker = yi, 
                                 subcoh = vi, 
                                 id = id, 
                                 data = tmpDat$data, 
                                 sets = tmpDat$Vi.k0.IND[,-c(1:2)], 
                                 predict.time = predict.time, 
                                 marker.cutpoint = cutoff, 
                                 estimation.method = ESTmethod, 
                                 Npert = 200)
    

    if(ESTmethod == "SP"){ 

      SE_ncc[i,] <- unlist(estimates_ncc$se)[c(1,2,4,3,6,5)] #change order to match up estimate names
    
      EST[i,] <- unlist(estimates_ncc$estimates)[c(1,2,4,3,6,5)]
      N_ncc[i] <-sum(tmpDat[[1]][,3]) # adding up across vi
      
    }else{
      
           
      SE_ncc[i,] <- unlist(estimates_ncc$se)[c(1,3,2,5,4)] #change order to match up estimate names
      
      EST[i,] <- unlist(estimates_ncc$estimates)[c(1,3,2,5,4)]
      N_ncc[i] <-sum(tmpDat[[1]][,3]) # adding up across vi
      
    }   
    
      
    if(i %%10==0) print(i)
  }
   
   

   SE_ncc <- as.data.frame(SE_ncc);
 
   EST <- as.data.frame(EST)

   
   names(SE_ncc) <- names(EST) <- mynames 

   out <- NULL
   
   for(tmppar in mynames ){

     if(tmppar =="AUC"){
       tmp.parval.H0 <- get.AUC.given.beta.intmethod(beta.H0, a=a, t = predict.time, f_x=f_x)
       tmp.parval.Ha <- get.AUC.given.beta.intmethod(beta.Ha, a=a, t = predict.time, f_x=f_x)
       
     }else if(tmppar=="beta"){
       
      tmp.parval.H0 <- beta.H0
      tmp.parval.Ha <- beta.Ha

      }else{
       fun.ind <- match(tmppar, c("TPF", "TPR", "FPR", "FPF", "NPV", "PPV"))
       
       myfun <- list(TPF.fun, TPF.fun, FPF.fun, FPF.fun, NPV.fun, PPV.fun)[[fun.ind]]
       
       tmp.parval.H0 <- myfun(c=cutoff, beta.H0, a=a, t = predict.time, f_x)
       tmp.parval.Ha <- myfun(c=cutoff, beta.Ha, a=a, t = predict.time, f_x)
     }
     
     tmpEST <- EST[,tmppar]
    
     tmpSE_ncc  <- SE_ncc[,tmppar]
         
   if(tmp.parval.H0 <= tmp.parval.Ha){
     
     power_ncc <- tmpEST - qnorm(1-alpha)*tmpSE_ncc > tmp.parval.H0

     logit.power_ncc <- logit(tmpEST) - qnorm(1-alpha)*((tmpSE_ncc)/(tmpEST*(1-tmpEST)))  > logit(tmp.parval.H0)

     out[[tmppar]] <- data.frame(cbind(est = tmpEST, 
                                       se_ncc = tmpSE_ncc, 
                                       power_ncc = power_ncc, 
                                       logit.power_ncc = logit.power_ncc))
     
   }else{
     
     power_ncc       <- tmpEST + qnorm(1-alpha)*tmpSE_ncc < tmp.parval.H0

     logit.power_ncc <- logit(tmpEST) + qnorm(1-alpha)*((tmpSE_ncc)/(tmpEST*(1-tmpEST)))  < logit(tmp.parval.H0)

     out[[tmppar]] <- data.frame(cbind(est = tmpEST, 
                                       se_ncc = tmpSE_ncc, 
                                       power_ncc = power_ncc, 
                                       logit.power_ncc = logit.power_ncc))
   }
   
   }
   #browser()
   if(ESTmethod == "SP"){
   out$beta$logit.power_ncc = out$beta$power_ncc
   }
  out$beta.H0 = beta.H0; out$beta.Ha = beta.Ha;
   
   out$a = a; 
  out$predict.time = predict.time; out$cutoff = cutoff; 
  out$N = N;out$N_ncc <- N_ncc;

   out$ESTmethod  <- ESTmethod; 
   out
}
