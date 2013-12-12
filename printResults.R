### function to print out the results of power simulations beautifully

#have three things
# 1. text
# 2. number
# 3. ggplot histogram
library(grid)
library(ggplot2)
#lets do it once and then build a function from there. 

# test 
system.time(
tst <- SimulateN( 1000,
          parameter= "TPR",
          S.0 = 0.9, 
          t.0 =1,
          cutoff=0, 
          predict.time = 2,
          parval.H0 = .5,
          parval.Ha =.55,
          f_x = dnorm, 
          F_xInv = qnorm,    
          mm = 100, 
          alpha = 0.05, 
          cens.perc = 20/100, 
          time.end = NULL )
)

tst

dev.off()
printResultPlot(tst, pars= c("beta", "AUC", "TPR", "FPR", "PPV", "NPV"), useLogit = FALSE)


## now work on how to display the results for the "model characteristics" page

# want to display:

#1. true values of all measures, including ranges 
#2. marker distribution along with cutpoint 
#3. Distribution of survival for Y at prediction time, and distribution of risk as a function of time at baseline
#4. ROC curve, and PPV/NPV curves

#again, do it once and then make a function from it
plotDistributions(
  parameter= "beta",
  S.0 = 0.9, 
  t.0 =1,
  cutoff=0 ,
  predict.time = 2,
  parval.H0 = .5,
  parval.Ha =.55)


grid.newpage()

### now I need a function to calculate ROC curve, PPV/NPV functions
tmp <- plotCurves(
  parameter= "AUC",
  S.0 = 0.9, 
  t.0 =1,

  predict.time = 2,
  parval.H0 = .5,
  parval.Ha =.55)


mytpf <- TPF.fun(0, .7, a = -log(.9), 2 ,f_x =dnorm)
myfpf <- FPF.fun(0,.7, a= -log(.9), 2 ,f_x =dnorm)
myppv <- PPV.fun(0, .7, a= -log(.9), 2 ,f_x =dnorm)
mynpv <- NPV.fun(0, .7, a= -log(.9), 2 ,f_x =dnorm)
pushViewport(viewport(layout = grid.layout(2,2, heights = unit(c(1,1), "null"))))

print(tmp$roc + geom_point(x=myfpf, y = mytpf, color = "blue", size = 5),
      vp = viewport(layout.pos.row =1:2, layout.pos.col = 1))

print(getCurvesPlot()$ppv+ geom_point(x=pnorm(cutoff), y = myppv, color = "blue", size = 5), 
      vp = viewport(layout.pos.row =1, layout.pos.col = 2))

print(getCurvesPlot()$npv+ geom_point(x=pnorm(cutoff), y = mynpv, color = "blue", size = 5), 
      vp = viewport(layout.pos.row =2, layout.pos.col = 2))


mydata <- SIM.data.singleMarker(nn = 500,
           mu = 0, 
           Sigma = 1, 
           beta = log(3), 
           lam0 = .1, cens.perc = .2, time.max = NULL)
 
names(mydata) = c("xi", "di", "Y")


getEstimates(data = mydata, 
                       cutoff,  
                       parameter = "PPV",
                       predict.time)





