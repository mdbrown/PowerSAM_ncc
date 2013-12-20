logit = function(x) log(x/(1-x))
expit = function(x) exp(x)/(1+exp(x))


############################3
## functions to plot results
#############################

printResultPlot <- function(simNobj, pars, useLogit = FALSE){
  # dev.off()
  index = 0;
  
  #open up the plot device
  pushViewport(viewport(layout = grid.layout(length(pars)+1,4, heights = unit(c(ifelse(length(pars)==1,3, 3),rep(10, length(pars))), "null"),
                                             widths = unit(c(1,1,1,2.5), "null")), gp = gpar(fontsize = 15)))
  
  grid.text("Measure \n (true value)", y = .6, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.text("Hypothesis Test",y = .5, x = 0, just = c("left", "center"), vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid.text("Power (%)",y = .5, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
  grid.text("Simulated Distribution", x = .4,  y = .5,vp=viewport(layout.pos.row = 1, layout.pos.col = 4))
  grid.lines(c(0,1), c(0.05,0.05), gp = gpar(col = "darkgrey"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:4))
  
  for(par in pars){
    index = index + 1; 
    
    #find the true values of the par
    if(par != "beta"){
      parTrue <- get.parameter.range(par, 
                                     simNobj$cutoff, 
                                     simNobj$a,
                                     simNobj$predict.time, f_x = dnorm, 
                                     lower = simNobj$beta.H0, 
                                     upper = simNobj$beta.Ha)
      
    }else{
      parTrue = round(c(simNobj$beta.H0, simNobj$beta.Ha), 2)
    }
    #prepare the data
    if(useLogit){
      x <- simNobj[[par]]
      x <- cbind(x$est, x$logit.power_ncc)
    }else{
      x <- simNobj[[par]]
      x <- cbind(x$est, x$power_ncc)
    }
    x<-as.data.frame(x)
    names(x) = c("estimate", "reject")
    
    #the ggplot histogram
    
    p <- ggplot(aes(estimate, fill = factor(reject)), data = x)
    p = p + geom_histogram(binwidth = diff(range(x$estimate, na.rm=TRUE))/30) + theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    ) 
    
    if(mean(x$reject, na.rm = TRUE)==1){
      #all reject 
      mycolors = c("#599ad3","#f9a65a")
    }else if(mean(x$reject, na.rm = TRUE)==0){
      #none reject
      mycolors = c("#f9a65a", "#599ad3")
    }else{
      mycolors <- c("#f9a65a", "#599ad3")
    }
    p = p + xlab("point estimate") + 
      scale_fill_manual(name = "test \nconclusion",
                        breaks = c("0", "1"), 
                        labels = c("fail to \nreject", "reject"),
                        values =  mycolors)
    p = p + geom_vline(xintercept = parTrue[2], size = 1.5, color ="#9e66ab")
    
    grid.text(parse(text = par), x = .5, y = .5, vp = viewport(layout.pos.row =index + 1, layout.pos.col = 1), gp = gpar(fontsize = 40))
    grid.text(paste("(=", round(parTrue[2], 2), ")"), x = .5, y = .3, vp =viewport(layout.pos.row = index+1, layout.pos.col = 1), gp = gpar(col = "#9e66ab", fontsize = 20, fontface = 2) )
    if(is.element(par, c("TPR", "FPR", "PPV", "NPV"))) par = paste(par, "(c)")
    
    tmp <- eval(substitute(paste("H[0]:",  mypar, myparsign, myparval),
                           list(mypar=par, myparsign =ifelse(par=="FPR", ">=", "<="), myparval = parTrue[1])))
    grid.text(parse(text = tmp), just = c("left", "center"), x = .05, y = .6, vp = viewport(layout.pos.row = index+1, layout.pos.col = 2))
    
    tmp <- eval(substitute(paste("H[a]:",  mypar, myparsign, myparval),
                           list(mypar=par, myparsign =ifelse(par=="FPR", "<", ">"), myparval = parTrue[1])))
    
    grid.text(parse(text = tmp), just = c("left", "center"), x = .05, y = .45, vp = viewport(layout.pos.row = index+1, layout.pos.col = 2))
    
    
    grid.text(paste(round(mean(x$reject, na.rm = TRUE), 2)*100), x = .5, y = .5,  gp = gpar(fontsize = 50, col = "#599ad3"), vp = viewport(layout.pos.row = index+1, layout.pos.col = 3))
    
    
    print(p, vp = viewport(layout.pos.row =index+1, layout.pos.col = 4))
    grid.lines(c(0,1), c(0.01,0.01), gp = gpar(col = "darkgrey"), vp = viewport(layout.pos.row = index+1, layout.pos.col = 1:4))
    
    
    
  }
  
}






plotDistributions <- function(
  parameter,
  S.0, 
  t.0 ,
  cutoff, 
  predict.time,
  parval.H0,
  parval.Ha, low.time=0, high.time=20){
  a = -log(S.0)/t.0
  
  
 

    if(parameter != "bet"){
      mybetas <- get.Betas(parameter, a, predict.time, cutoff, parval.H0, parval.Ha)
      
      parTrue <- round(c(parval.H0, parval.Ha), 2)
    }else{
      mybetas <-c(parval.H0, parval.Ha)
      parTrue = round(c(parval.H0, parval.Ha), 2)
    }
  
  ##lets do marker distribution with cutpoint along wit distribution of risk 
  
  
  Yvals <- seq(-3, 3, by = .05)
  MarkerDist <- data.frame(cbind(Y = Yvals, fy = dnorm(Yvals), risk = Yvals > cutoff))
  .e <- environment()
  Ydist_p <- ggplot(data = MarkerDist, aes(x=Y, y = fy), environment = .e) +
    geom_area(aes(fill = factor(risk))) + 
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    )  + xlab("Marker Y") + ylab("f(Y)")+ theme(legend.position = "bottom") + 
    geom_vline(xintercept = cutoff, col ="#9e66ab", size = 1.5 )+
    geom_text(aes(ifelse(cutoff>0, -2, 2), .3,
                  label = paste("Marker \n cutpoint = ", cutoff)) , color ="#9e66ab", size = 5)+
    scale_fill_manual(name = "",
                      breaks = c("0", "1"), 
                      labels = c("low risk", "high risk"),
                      values =  c("#599ad3","#f9a65a" )) +
    geom_text(aes(-2.5, .1, label = paste(round(pnorm(cutoff), 2)*100, "%", sep ="" )) , color ="#599ad3", size = 6) +
    geom_text(aes( 2.5, .1, label = paste(round(1-pnorm(cutoff), 2)*100, "%", sep ="" )) , color ="#f9a65a", size = 6)
  
  
  
  
  ##now for distribution of S(t |Y)
  

  Yvals <- round(seq(-2.5, 2.5, length.out = 9), 1)
  tvals <- seq(low.time, high.time, length.out = 25)


  SyDist <- NULL
  for(Y in Yvals){
    SyDist = c(SyDist, Sy_fun(-log(S.0)/t.0, Y, tvals, mybetas[2]))
    
  }
  
  
  SyDist = data.frame(Sy = SyDist, t = rep(tvals, length(Yvals)), Y = rep(Yvals, rep(length(tvals), length(Yvals))))

  
  Sy_p <- ggplot(aes(x = t, y = Sy, color = factor(Y)), data = SyDist) +
    geom_line(size = 2) + geom_vline(xintercept = predict.time, color = "#525252", size= 1.8, linetype = 2) +
   #   geom_text(aes(30, .75, label = as.character(expression(S(t,Y) == S[0](t)^{e^{-beta*Y}}))), parse = TRUE, color = "black", size = 8)+
    scale_color_manual(name = "Value \nof Y", values = brewer.pal(9, "Spectral")) +
    xlab("Time t") +ylab("S(t |Y) = Pr(T > t | Y)")
  
  
  #now Sy by Y at predict.time
  
  Yvals = seq(-3, 3, by = .1)
  
  SyDist = c( Sy_fun(a = -log(S.0)/t.0, Y= Yvals, predict.time, mybetas[2]))
  SyDist = data.frame(Sy = SyDist, Y = Yvals)  
  
  SybyY_p <- ggplot(aes(x = Y, y = Sy), data = SyDist) + 
    geom_line(size = 1.8, colour ="#525252" , linetype = 2) +theme(panel.border = element_blank()) + xlab("Marker Y") + 
    geom_vline(xintercept = cutoff, col ="#9e66ab", size = 1.5 ) +
    ylab(paste("P( T >",predict.time, "| Y)") ) +ylim(0,1)

  pushViewport(viewport(layout = grid.layout(2,2, heights = unit(c(2,3), "null"))))
  
  print(Ydist_p, vp = viewport(layout.pos.row =2, layout.pos.col = 2))
  
  print(Sy_p, vp = viewport(layout.pos.row =c(1:2), layout.pos.col = 1))
  print(SybyY_p, vp = viewport(layout.pos.row =1, layout.pos.col = 2))
  
}

Sy_fun <- function(a, Y, t, beta){
  exp(-a*t*exp(Y*beta))
}



#first ROC
plotCurves <- function( parameter, 
  S.0 , 
  t.0 ,
  cutoff, 
  predict.time ,
  parval.H0, 
  parval.Ha
){
  a = -log(S.0)/t.0
  if(parameter != "beta"){
    mybetas <- get.Betas(substr(parameter, 1,3), a, predict.time, cutoff = cutoff, parval.H0, parval.Ha)
  }else{
    mybetas <-c(parval.H0, parval.Ha)
  }
  
  
  beta.H0 <- mybetas[1]
  beta.Ha <- mybetas[2]
  c.points <- qnorm(1:99/100)
  tpf <- numeric(length(c.points));
  fpf <- numeric(length(c.points)); 
  ppv <- numeric(length(c.points))
  npv <- numeric(length(c.points))
  a <- -log(S.0)/t.0
  i=0
  for(c in c.points){
    i = i+1
    tpf[i] <- TPF.fun(c, beta.Ha, a, predict.time ,f_x =dnorm)
    fpf[i] <- FPF.fun(c, beta.Ha, a, predict.time ,f_x =dnorm)
    ppv[i] <- PPV.fun(c, beta.Ha, a, predict.time ,f_x =dnorm)
    npv[i] <- NPV.fun(c, beta.Ha, a, predict.time ,f_x =dnorm)
    
  }

 # mytpf <- TPF.fun(cutoff, beta.Ha, a, predict.time ,f_x =dnorm)
#  myfpf <- FPF.fun(cutoff, beta.Ha, a, predict.time ,f_x =dnorm)
#  myppv <- PPV.fun(cutoff, beta.Ha, a, predict.time ,f_x =dnorm)
#  mynpv <- NPV.fun(cutoff, beta.Ha, a, predict.time ,f_x =dnorm)
  
  roc.dat.Ha <- data.frame(c = c.points, tpf = tpf, fpf = fpf, ppv = ppv, npv = npv)
  
  
  i=0
  for(c in c.points){
    i = i+1
    tpf[i] <- TPF.fun(c, beta.H0, a, predict.time ,f_x =dnorm)
    fpf[i] <- FPF.fun(c, beta.H0, a, predict.time ,f_x =dnorm)
    ppv[i] <- PPV.fun(c, beta.H0, a, predict.time ,f_x =dnorm)
    npv[i] <- NPV.fun(c, beta.H0, a, predict.time ,f_x =dnorm)
    
  }
  
  roc.dat.H0 <- data.frame(c = c.points, tpf = tpf, fpf = fpf, ppv = ppv, npv = npv)
  roc.dat <- rbind(roc.dat.H0, roc.dat.Ha)
  roc.dat$Hypothesis <- factor(c(rep("H0", 99),rep("Ha", 99)))
  roc.dat$quantile.c <- pnorm(c.points)
  roc <- ggplot(data = roc.dat, aes(x = fpf, y = tpf, colour = Hypothesis)) 
  roc <- roc+geom_line(size = 1.5)+
    #geom_abline(a=0, b = 1, linetype = 2, size = 1.5) +
    #geom_point(x=myfpf, y = mytpf, color = "blue", size = 5) + 
    xlim(0,1) + ylim(0,1) + scale_colour_discrete("", breaks = c("H0", "Ha"), labels =c("Value under \n    Null", "True Value")) +
    ylab("TPR(c)") + xlab( "FPR(c)")+
    theme(legend.position = "top", legend.text = element_text(size = 13))
  
  ppv <-  ggplot(data = roc.dat, aes(x = quantile.c, y = ppv, colour = Hypothesis)) +
    geom_line(size = 1.5) + 
    #geom_point(x=pnorm(cutoff), y = myppv, color = "blue", size = 5) + 
    xlim(0,1) + ylim(0,1)+ xlab("F(c)") + ylab("PPV(c)") +
   scale_colour_discrete("", breaks = c("H0", "Ha"), labels =c("Value under \n    Null", "True Value")) +
    theme(legend.position = "none")
  
  npv <-  ggplot(data = roc.dat, aes(x = quantile.c, y = npv, colour = Hypothesis)) +
    geom_line(size = 1.5) + 
    #geom_point(x=pnorm(cutoff), y = mynpv, color = "blue", size = 5) + 
    xlim(0,1) + ylim(0,1) + xlab("F(c)") + ylab("NPV(c)")+
    scale_colour_discrete("", breaks = c("H0", "Ha"), labels =c("Value under \n    Null", "True Value"))+
    theme(legend.position = "none")
  
  
  


  
  mytpf <- TPF.fun(cutoff, mybetas[2], a, predict.time ,f_x =dnorm)
  myfpf <- FPF.fun(cutoff, mybetas[2], a, predict.time ,f_x =dnorm)
  myppv <- PPV.fun(cutoff, mybetas[2], a, predict.time ,f_x =dnorm)
  mynpv <- NPV.fun(cutoff, mybetas[2], a, predict.time ,f_x =dnorm)
  
  nutpf <- TPF.fun(cutoff, mybetas[1], a, predict.time ,f_x =dnorm)
  nufpf <- FPF.fun(cutoff, mybetas[1], a, predict.time ,f_x =dnorm)
  nuppv <- PPV.fun(cutoff, mybetas[1], a, predict.time ,f_x =dnorm)
  nunpv <- NPV.fun(cutoff, mybetas[1], a, predict.time ,f_x =dnorm)
  
  
  pushViewport(viewport(layout = grid.layout(3,1, heights = unit(c(2,1, 1), "null"))))
  
  print(roc +
          geom_point(x=myfpf, y = mytpf, color = "#00BFC4", size = 4) +
          geom_point(x=nufpf, y = nutpf, color = "#F8766D", size = 4),
        vp = viewport(layout.pos.row =1, layout.pos.col = 1))
  
  print(ppv+ 
          geom_point(x=pnorm(cutoff), y = myppv, color = "#00BFC4", size = 4)+
          geom_point(x=pnorm(cutoff), y = nuppv, color = "#F8766D", size = 4), 
        vp = viewport(layout.pos.row =2, layout.pos.col = 1))
  
  print(npv+
          geom_point(x=pnorm(cutoff), y = mynpv, color = "#00BFC4", size = 4)+
          geom_point(x=pnorm(cutoff), y = nunpv, color = "#F8766D", size = 4), 
        vp = viewport(layout.pos.row =3, layout.pos.col = 1))
  
  
}




sum.I<-function(yy,FUN,Yi,Vi=NULL){
  
  if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}
  
  pos <- rank(c(yy,Yi),ties.method='f')[1:length(yy)]-rank(yy,ties.method='f')
  
  if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos
  
  if (!is.null(Vi)) {
    
    if(substring(FUN,2,2)=="=") tmpind <- order(-Yi) else  tmpind <- order(Yi)
    
    Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
    
    return(rbind(0,Vi)[pos+1,])
    
  } else return(pos)
}

##################################
## Calculate true measure values and range
##################################

get.Betas <- function(parameter, a, predict.time, cutoff, parval.H0, parval.Ha,  f_x = dnorm ){
if(parameter =="AUC"){

  beta.H0 <- uniroot(function(beta, a, t, f_x, AUC) get.AUC.given.beta.intmethod(beta, a, t, f_x=f_x) - AUC, 
                     lower = -10, upper = 100,  a = a, t = predict.time,
                     f_x=f_x, AUC = parval.H0)$root
  
  
  beta.Ha <- uniroot(function(beta, a, t, f_x, AUC) get.AUC.given.beta.intmethod(beta, a, t, f_x=f_x) - AUC, 
                     lower = -10, upper = 100,  a = a, t = predict.time,
                     f_x=f_x, AUC = parval.Ha)$root
  
}else{


  fun.ind <- match(parameter, c("TPF", "TPR", "FPR", "FPF", "NPV", "PPV"))
  
  myfun <- list(TPF.fun, TPF.fun, FPF.fun, FPF.fun, NPV.fun, PPV.fun)[[fun.ind]]
  
  #find beta.H0/a: the beta that makes parval.H0/a
  beta.H0 <- uniroot(function(beta, c, a, t, f_x, parval) myfun(c, beta, a, t, f_x) - parval, 
                     lower = -10, upper = 100,   
                     c = cutoff,  a = a, t = predict.time,
                     f_x=f_x, parval = parval.H0)$root
  
  beta.Ha <- uniroot(function(beta, c, a, t, f_x, parval) myfun(c, beta, a, t, f_x) - parval, 
                     lower = -10, upper = 100,   
                     c = cutoff,  a = a, t = predict.time,
                     f_x=f_x, parval = parval.Ha)$root
  
  
  
  
}
c(beta.H0, beta.Ha)
}




get.parameter.range <- function(parameter, cutoff, a, predict.time, f_x, lower = 0, upper = 75){

  if(parameter != "AUC"){

      fun.ind <- match(parameter, c("TPF", "TPR", "FPR", "FPF", "NPV", "PPV"))

      myfun <- list(TPF.fun, TPF.fun, FPF.fun, FPF.fun, NPV.fun, PPV.fun)[[fun.ind]]
     # browser()
      #find beta.Ha: the beta that makes parval.H0
      upper.par <- myfun(c = cutoff, beta=upper, a = a, t=predict.time, f_x=f_x)   
      lower.par <- myfun(c = cutoff, beta=lower, a = a, t=predict.time, f_x=f_x)

  }else{

      upper.par <- get.AUC.given.beta.intmethod(beta = upper, a = a, t = predict.time, f_x=f_x) 
      lower.par <- get.AUC.given.beta.intmethod(beta = lower, a = a, t = predict.time, f_x=f_x) 

  }
 if(lower.par <0) lower.par = 0.00
 if(upper.par>=.99) upper.par = .99

round(c(lower.par, upper.par), 4) 

}




## here we are explicitely calculating the double integral to get AUC
## this works fastest of the methods
get.AUC.given.beta.intmethod <- function( beta, a, t, f_x ){

#get Pr(No Disease) = Pr(T>t) = integral (s(t|x))*f_x
int.x <- function(x, beta = beta,  a = a, t = t, f_x = f_x)  (exp(-a*t*exp(beta*x)))*f_x(x)
Pr.noD  <- integrate( int.x, lower = -Inf, upper =Inf, beta = beta,  a = a, t = t, f_x = f_x)$value

Pr.D <- 1-Pr.noD
  

myfun <- function(y, x, beta, a, t, f_x) (1-exp(-a*t*exp(beta*x)))*f_x(x)*(exp(-a*t*exp(beta*y)))*f_x(y)
#double integral to calculate Pr(x_disease<x_nodisease) 
out <- integrate(function(x) {
         sapply(x, function(mu) {
             
             integrate(myfun, lower = -Inf, upper = mu,  x = mu, beta = beta, a = a, t = t, f_x =f_x, abs.tol = 1e-8)$value
             })
        }, -Inf, Inf)$value

out <- out/(Pr.noD*Pr.D)
#Pr(x_disease>x_nodisease) = AUC 
out


}

###

TPF.fun <- function(c, beta, a, t, f_x){

  #function to integrate from c to infinity over (for the numerator)
  # the denominator is the integral over the whole domain (so it is equal to Pr(T<=t) = 1-S(t)
  # this is Pr(T<=t, | Y = y)*Pr(Y=y) = (1-S( t| y))*f(y) 
  

  int.x = function(x, beta = beta,  a = a, t = t, f_x = f_x)  (1-exp(-a*t*exp(beta*x)))*f_x(x)
  
  num <- integrate( int.x, lower = c   , upper =Inf, beta = beta,  a = a, t = t, f_x = f_x)$value
  den <- integrate( int.x, lower = -Inf, upper =Inf, beta = beta,  a = a, t = t, f_x = f_x)$value
  num/den
}



FPF.fun <- function(c, beta, a, t ,f_x){

  #function to integrate from c to infinity over (for the numerator)
  # this is Pr(T>t, | Y = y)*Pr(Y=y) = (S( t| y))*f(y) 
  int.x = function(x, beta = beta,  a = a, t = t, f_x = f_x)  (exp(-a*t*exp(beta*x)))*f_x(x)
  
  num <- integrate( int.x, lower = c  , upper =Inf, beta = beta, a = a, t = t, f_x = f_x)$value
  den <- integrate( int.x, lower =-Inf, upper =Inf, beta = beta, a = a, t = t, f_x = f_x)$value
  
  num/den
}



NPV.fun <- function(c, beta, a, t, f_x){

  # this is Pr(T>t, | Y = y)*Pr(Y=y) = (S( t| y))*f(y) 
  int.x <- function(x, beta = beta,  a = a, t = t, f_x = f_x)  (exp(-a*t*exp(beta*x)))*f_x(x)

  num <- integrate( int.x, lower = -Inf, upper = c, beta = beta, a = a, t = t, f_x = f_x)$value
  den <- integrate( f_x,   lower = -Inf, upper = c)$value

  num/den

}


PPV.fun <- function(c, beta, a, t, f_x){
  # this is Pr(T>t, | Y = y)*Pr(Y=y) = (S( t| y))*f(y) 
  int.x = function(x, beta = beta,  a = a, t = t, f_x = f_x)  (1-exp(-a*t*exp(beta*x)))*f_x(x)
  
  num <- integrate( int.x, lower = c, upper = Inf, beta = beta, a = a, t = t, f_x = f_x)$value
  den <- integrate( f_x,   lower = c, upper = Inf)$value


  num/den

}





