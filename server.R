require(shiny)
require(survival)
require(ggplot2)
require(shinyIncubator)



require(grid)
require(RColorBrewer)
#options(markdown.HTML.options=c('fragment_only','base64_images'))
source("global.R")
source("main.R")
source("subroutines.R")
source("Estimation.R")
source("CoxPtb.R")

# Define server logic required to summarize and view the selected dataset
shinyServer(function(input, output) {

  
  ### contents 
  
  # 1. Reactive Functions
  # 2. render UI functions
  # 3. render text
  # 4. render tables
  # 5. render plots
  
  
#####
## Reactive Functions
#####
  
   # simulate power for a given sample size
 getSimulateN <- reactive({
    if(input$runSim!=0){

     isolate({ SimulateN(N = input$cohortN,
             parameter= substr(input$parameter, 1, 3),
             S.0 = input$S.0, 
             t.0 = input$t.0,
             cutoff=input$cutoff, 
             predict.time = input$predict.time,
             parval.H0 = input$H0,
             parval.Ha = input$Ha,
             ESTmethod = input$ESTmethod, 
             f_x = dnorm, 
             F_xInv = qnorm,    
             mm = input$mmSim,  
             alpha = input$alpha, 
             cens.perc = input$cens.perc/100, 
             time.end = NULL, time.max = input$time.max, censorType= input$censorType,
             type = input$typeOfStudy, 
             nmatch = input$nmatch)
     })
    
    }
  })

# get parameter ranges for different risk measures
  getParRange <- reactive({

   get.parameter.range(parameter = substr(input$parameter, 1,3),
                       cutoff    = input$cutoff, 
                       a         = -log(input$S.0)/input$t.0, 
                       predict.time   = input$predict.time, 
                       f_x       = dnorm)    


})
 
  

PrNocens <- reactive({
  if(input$censorType=="cens.perc") Pr.Nocens = 1-input$cens.perc/100
  else {
    a = -log(input$S.0)/input$t.0
    
    my_fun <- function( Y,a, t, beta){
      (1- exp(-a*t*exp(Y*beta)))*dnorm(Y)
    }
    if(substr(input$parameter, 1,3)!="bet"){
      mybetas <- get.Betas(substr(input$parameter, 1,3), a, input$predict.time, cutoff = input$cutoff, input$H0, input$Ha)
    }else{
      mybetas <-c(input$H0, input$Ha)
    }
    
    Pr.Nocens <- integrate(my_fun, lower = -Inf, upper = Inf, a= a, t= input$time.max, beta = mybetas[2] )$value
    
  }
  
  Pr.Nocens
  
})


#####
## render ui's
#####
 
  
  output$censType <- renderUI({
    if(input$censorType=="cens.perc") sliderInput("cens.perc", label = "Percentage Censored:", min = 0, max = 99, value = 50 )
    else numericInput("time.max", label = "", value=10, min = 0)
  })
  
  output$NullLine1 <- reactive({ 
    HTML(paste("H<sub>0</sub>: ", input$parameter,  ifelse(input$parameter=="FPR(c)", " &ge; ", " &le; "), input$H0, " vs.", sep = ""))
  }) 
  output$NullLine2 <- reactive({
    HTML(paste("H<sub>a</sub>: ", input$parameter,  ifelse(input$parameter=="FPR(c)", " < ", " > "),input$H0, sep = ""))
  })
  
  output$nullInput <- renderUI({
    
    sliderInput("H0", label = "", 
                min = ifelse(input$parameter =="FPR(c)", getParRange()[2], getParRange()[1]), 
                max = ifelse(input$parameter =="FPR(c)", getParRange()[1], getParRange()[2]),
                value = ifelse(input$parameter =="FPR(c)", max(getParRange()), min(getParRange())),
                step = .005)
    
  })
  
  
  
  output$altInput <- renderUI({
    
    sliderInput("Ha", label = paste("expected value of", input$parameter),
                min = ifelse(input$parameter =="FPR(c)", getParRange()[2], input$H0), 
                max = ifelse(input$parameter =="FPR(c)", input$H0, getParRange()[2]),
                value = mean(input$H0, getParRange()[2]),
                step = .005)
    
  })


#####
## Render text 
#####



 output$printEstMethod <-renderText({
   
   if(is.null(getSimulateN())) return(NULL)
   
   
 if(getSimulateN()$ESTmethod=="SP") ESTmethod = "Semi-parametric"
 else ESTmethod = "Non-parametric"
   HTML(paste("Estimation Method: <strong>", ESTmethod, "</strong>", sep = ""))
   
 }) 
 
  output$printSampSize <-renderText({

      if(is.null(getSimulateN())) return(NULL)

    HTML(paste("Cohort sample size: <strong>", getSimulateN()$N, "</strong>", sep = ""))
    
  }) 
 
 output$printSampSizeNCC <-renderText({
   if(is.null(getSimulateN())) return(NULL)
   HTML(paste(ifelse(getSimulateN()$type ==1, "(Average) ", " "), "NCC sample size: <strong>", round(mean(getSimulateN()$N_ncc)), "</strong>" ,sep = ""))
   
 }) 
  

output$getEventRates <- renderText({
  a = -log(input$S.0)/input$t.0
  
  my_fun <- function( Y,a, t, beta){
    (1- exp(-a*t*exp(Y*beta)))*dnorm(Y)
  }
  if(substr(input$parameter, 1,3)!="bet"){
    mybetas <- get.Betas(substr(input$parameter, 1,3), a, input$predict.time, cutoff = input$cutoff, input$H0, input$Ha)
  }else{
    mybetas <-c(input$H0, input$Ha)
  }
  
  Event.rate <- integrate(my_fun, lower = -Inf, upper = Inf, a= a, t= input$predict.time, beta = mybetas[2] )$value
  
  if(input$censorType=="cens.perc") {
    Pr.Nocens  = 1-input$cens.perc/100
    out = Event.rate*Pr.Nocens 
  }
  else{
    
    if(input$predict.time <= input$time.max) out = Event.rate
    else out =integrate(my_fun, lower = -Inf, upper = Inf, a= a, t= input$time.max, beta = mybetas[2] )$value
    
  }
  
  x <- HTML(paste("<li><strong>Event rate at prediction time: Pr( T < prediction time & T not censored) =",  round(out, 3), "</strong></li>"))
  x
  
})

output$censoringPercentage <- renderText({
  
  
  
  if(input$censorType=="cens.perc") Pr.Nocens = 1-input$cens.perc/100
  else {
    a = -log(input$S.0)/input$t.0
    
    my_fun <- function( Y,a, t, beta){
      (1- exp(-a*t*exp(Y*beta)))*dnorm(Y)
    }
    if(substr(input$parameter, 1,3)!="bet"){
      mybetas <- get.Betas(substr(input$parameter, 1,3), a, input$predict.time, cutoff = input$cutoff, input$H0, input$Ha)
    }else{
      mybetas <-c(input$H0, input$Ha)
    }
    
    Pr.Nocens <- integrate(my_fun, lower = -Inf, upper = Inf, a= a, t= input$time.max, beta = mybetas[2] )$value
    
  }
  
  x <- HTML(paste("<li><strong>Overall observed event rate (due to censoring): Pr( T not censored) =", round(Pr.Nocens,3), "</strong></li>"))
  x
  
})


output$NCCsamplesizeALL <- renderText({

  a = -log(input$S.0)/input$t.0
  
  if(substr(input$parameter, 1,3)!="bet"){
    mybetas <- get.Betas(substr(input$parameter, 1,3), a, input$predict.time, cutoff = input$cutoff, input$H0, input$Ha)
  }else{
    mybetas <-c(input$H0, input$Ha)
  }
  
  
  browser()
  tmpDat <- SIM.data.singleMarker(nn = input$cohortN, 
                                  beta = mybetas[2], 
                                  lam0 = a, 
                                  cens.perc = input$cens.perc/100, 
                                  time.max = input$time.max, 
                                  m.match = input$nmatch)   
  
  
  subdata <- data.frame(tmpDat[[1]])
  subdata <- subdata[subdata$vi==1,] 
  
  weights <- FNCC.WGT.FUN(tmpDat[[1]],
                          V.IND = tmpDat[[2]], 
                          Iik0 = tmpDat[[3]], 
                          Zmatch.ind=NULL, 
                          m.match = input$nmatch)
    
  trueW <- weights[subdata$di==0] 
  
  tmp1 <- with(subdata, sapply(xi[di==0], function(x) sum(x < xi[di==1])))
  tmp2 <- with(subdata, sapply(xi[di==0], function(x) sum(x < xi)))
  trueW[1:10]
  (tmp1/tmp2)[1:10]
  
  tmpN = min( input$cohortN*(Pr.Nocens) + input$cohortN*(1-Pr.Nocens)*input$nmatch*tmp.prob, 
              input$cohortN)
  
  
  
  
  HTML(paste("<h3>", "Average total sample size under NCC design: <strong>", tmpN, "</strong></h3>", sep = ""))
  
})


######
## Render tables
######



output$NCCsamplesize <- renderTable({
  
  out = data.frame(cbind(c(1,2), c(3,4)))
  
  
  
  
  names(out) <- c("NCC sample", "Cohort")
  row.names(out) <- c("Controls", "Cases" )
  
  round(out)
}, digits = 0)


output$TrueValuesTable <- renderTable({
  
  a = -log(input$S.0)/input$t.0
  
  if(substr(input$parameter, 1,3)!="bet"){
    mybetas <- get.Betas(substr(input$parameter, 1,3), a, input$predict.time, cutoff = input$cutoff, input$H0, input$Ha)
  }else{
    mybetas <-c(input$H0, input$Ha)
  }
  
  
  altauc <- get.AUC.given.beta.intmethod(mybetas[2], a, input$predict.time, f_x=dnorm )
  alttpf <- TPF.fun(input$cutoff, mybetas[2], a, input$predict.time ,f_x =dnorm)
  altfpf <- FPF.fun(input$cutoff, mybetas[2], a, input$predict.time ,f_x =dnorm)
  altppv <- PPV.fun(input$cutoff, mybetas[2], a, input$predict.time ,f_x =dnorm)
  altnpv <- NPV.fun(input$cutoff, mybetas[2], a, input$predict.time ,f_x =dnorm)
  
  nullauc <- get.AUC.given.beta.intmethod(mybetas[1], a, input$predict.time, f_x=dnorm )
  nulltpf <- TPF.fun(input$cutoff, mybetas[1], a, input$predict.time ,f_x =dnorm)
  nullfpf <- FPF.fun(input$cutoff, mybetas[1], a, input$predict.time ,f_x =dnorm)
  nullppv <- PPV.fun(input$cutoff, mybetas[1], a, input$predict.time ,f_x =dnorm)
  nullnpv <- NPV.fun(input$cutoff, mybetas[1], a, input$predict.time ,f_x =dnorm)
  
  mytable <- data.frame(rbind( c(round(abs(mybetas), 3)), 
                               c(nullauc, altauc), 
                               c(nulltpf, alttpf), 
                               c(nullfpf, altfpf), 
                               c(nullppv, altppv), 
                               c(nullnpv, altnpv)))
  names(mytable ) = c("Value under Null Hypothesis", "True Value")
  
  row.names(mytable) = c("β", "AUC", 
                         paste("TPR(", input$cutoff ,")", sep = ""), 
                         paste("FPR(", input$cutoff ,")", sep = ""),
                         paste("PPV(", input$cutoff ,")", sep = ""), 
                         paste("NPV(", input$cutoff ,")", sep = ""))
  mytable
})



#### 
## Render plots
####
 
  output$DistributionPlot <- renderPlot({
    if(input$plotNow!=0){
     return(isolate(plotDistributions(parameter = substr(input$parameter, 1,3),
                      S.0 = input$S.0, 
                      t.0 = input$t.0 ,
                      cutoff = input$cutoff, 
                      predict.time = input$predict.time,
                      parval.H0 = input$H0,
                      parval.Ha = input$Ha,
                       low.time = 0,#input$lowertime, 
                       high.time = input$uppertime)))
    }

    })
  

  
  output$CurvesPlot <- renderPlot({
      if(input$plotNow!=0){
        
        return(isolate(plotCurves(parameter = substr(input$parameter, 1,3),
                                  S.0 = input$S.0, 
                                  t.0 = input$t.0 ,
                                  cutoff = input$cutoff, 
                                  predict.time = input$predict.time,
                                  parval.H0 = input$H0, 
                                  parval.Ha = input$Ha)))
      }
  })
  
  
 


#main plot for Simulation
output$simulationGraphTop <- renderPlot({
  if(!is.null(getSimulateN())){
    
    printResultPlot(getSimulateN(), pars= substr(input$parameter, 1,3), useLogit = FALSE)
  }
})

output$simulationGraphSub <- renderPlot({
  if(!is.null(getSimulateN())){
    
    if(getSimulateN()$ESTmethod=="SP") mypars = c( "beta", "AUC", "TPR", "FPR", "PPV", "NPV")
    else mypars = c( "AUC", "TPR", "FPR", "PPV", "NPV")
    
    
    mypars = mypars[mypars!=substr(input$parameter, 1,3)]
    if(length(mypars)>0){
      printResultPlot(getSimulateN(), pars= mypars, useLogit = FALSE)
    }else{ return()}
    
  }
})

  
  
  
})
