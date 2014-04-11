require(shiny)
require(survival)
require(ggplot2)
require(shinyIncubator)
require(survMarkerTwoPhase)


require(grid)
require(RColorBrewer)
#options(markdown.HTML.options=c('fragment_only','base64_images'))

source("SimulationFuns.R")
source("subroutines.R")


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
  
## 
getBetas <- reactive({
  a = -log(input$S.0)/input$t.0
  
  
  if(substr(input$parameter, 1,3)!="bet"){
    mybetas <- get.Betas(substr(input$parameter, 1,3), a, input$predict.time, cutoff = input$cutoff, input$H0, input$Ha)
  }else{
    mybetas <-c(input$H0, input$Ha)
  }
  
  mybetas
  
})
getSampleSizes <- reactive({
      
      a = -log(input$S.0)/input$t.0
      
      if(substr(input$parameter, 1,3)!="bet"){
        mybetas <- get.Betas(substr(input$parameter, 1,3), 
                             a, 
                             input$predict.time, 
                             cutoff = input$cutoff, 
                             input$H0, input$Ha)
      }else{
        mybetas <-c(input$H0, input$Ha)
      }
      
      Lam = getLam()
      tmpDat <- data.frame(SIM.data.singleMarker(nn = input$cohortN, 
                                      beta = mybetas[2], 
                                      lam0 = a, 
                                      cens.lam = Lam, 
                                      time.max = input$time.max, 
                                      m.match = input$nmatch)[[1]])

      out <- with(tmpDat, table(vi, di))
      out
})


   # simulate power for a given sample size
 getPowerSim <- reactive({
    if(input$runSim!=0){
      myLam = getLam()
     isolate({ PowerSim(N = input$cohortN,
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
             cens.lam = myLam, 
             time.max = input$time.max,
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
 
  

EventRates <- reactive({
  
  a = -log(input$S.0)/input$t.0
  
  beta = getBetas()[2]
  lam = getLam()
  
  
  eventrate <- Ft(input$predict.time, a = a, beta) - (Pr.cTtmax(lam, tmax =input$predict.time , a, beta))
  censoringrate <-  Ft(input$time.max, a = a, beta) - (Pr.cTtmax(lam, tmax =input$time.max , a, beta))
  list("eventrate" = eventrate, "censoringrate" = censoringrate)
  
})




output$NullLine1 <- reactive({ 
  HTML(paste("H<sub>0</sub>: ", input$parameter,  ifelse(input$parameter=="FPR(c)", " &ge; ", " &le; "), input$H0, " vs.", sep = ""))
}) 
output$NullLine2 <- reactive({
  HTML(paste("H<sub>a</sub>: ", input$parameter,  ifelse(input$parameter=="FPR(c)", " < ", " > "),input$H0, sep = ""))
})


getLam <- reactive({
  a = -log(input$S.0)/input$t.0
  beta = getBetas()[2]
  mymax = 100-floor(Pr.cTtmax(lam = 100, input$time.max, a, getBetas()[2])*100)
  if(input$cens.perc >mymax){
    
    
    out <- uniroot(function(lam, cens.perc, tmax, a, beta){Pr.cTtmax(lam, tmax, a, beta) - cens.perc}, 
                   c(0, 100),
                   cens.perc =  (input$cens.perc - mymax)/100, 
                   tmax = input$time.max, 
                   a = a, 
                   beta = beta)$root
  }else{
    
    out = 0
  }
  
  out
  
})

#####
## render ui's
#####
 

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
output$censoringInput <- renderUI({
  
  a = -log(input$S.0)/input$t.0
  mymax = 100-floor(Pr.cTtmax(lam = 100, input$time.max, a, getBetas()[2])*100)
  #mymax = floor((1- Ft(input$time.max, a = a, getBetas()[[2]]))*100)
  
  sliderInput("cens.perc",
              "Total percent of observations censored:", 
              min = mymax, 
              max = 99, 
              value =mymax , 
              step=1)
  
  
})


#####
## Render text 
#####


output$censoringNote <- renderText({
  a = -log(input$S.0)/input$t.0
  paste("This followup time would cause ", 
        100-floor(Pr.cTtmax(lam = 100, input$time.max, a, getBetas()[2])*100), 
        "% of observations to be censored on average. To add censoring due to study participants dropping out before time to followup, use the slider below.", sep = "")
})

 output$printEstMethod <-renderText({
   
   if(is.null(getPowerSim())) return(NULL)
   
   
 if(getPowerSim()$ESTmethod=="SP") ESTmethod = "Semi-parametric"
 else ESTmethod = "Non-parametric"
   HTML(paste("Estimation Method: <strong>", ESTmethod, "</strong>", sep = ""))
   
 }) 
 
  output$printSampSize <-renderText({

      if(is.null(getPowerSim())) return(NULL)

    HTML(paste("Cohort sample size: <strong>", getPowerSim()$N, "</strong>", sep = ""))
    
  }) 
 
 output$printSampSizeNCC <-renderText({
   if(is.null(getPowerSim())) return(NULL)
   HTML(paste(ifelse(getPowerSim()$type ==1, "(Average) ", " "), "NCC sample size: <strong>", round(mean(getPowerSim()$N_ncc)), "</strong>" ,sep = ""))
   
 }) 
  

output$getEventRates <- renderText({
  
  x <- HTML(paste("<li><strong>Event rate at prediction time: Pr( T < prediction time & T not censored) =",  round(EventRates()$eventrate, 2), "</strong></li>"))
  x
  
  
})

output$censoringPercentage <- renderText({
  
  x <- HTML(paste("<li><strong>Overall observed event rate (due to censoring): Pr( T not censored) =", round(EventRates()$censoringrate,2), "</strong></li>"))
  x
  
})


output$NCCsamplesizeALL <- renderText({
  

  HTML(paste("<p>", "A simulated data set with the input parameters provided had NCC sample size of: <strong>", sum(getSampleSizes()[2,]), " observations</strong></p>", sep = ""))
  
})


######
## Render tables
######



output$NCCsamplesize <- renderTable({
  
  out = matrix(0, nrow = 2, ncol = 2)
  out[,1] = getSampleSizes()[2,]
  out[,2] = getSampleSizes()[1,] + getSampleSizes()[2,]
  
  out <- as.data.frame(out)
  
  
  
  names(out) <- c("NCC sample", "Cohort")
  row.names(out) <- c("Controls", "Cases" )
  
  out
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
  if(!is.null(getPowerSim())){
    
    printResultPlot(getPowerSim(), pars= substr(input$parameter, 1,3), useLogit = FALSE)
  }
})

output$simulationGraphSub <- renderPlot({
  if(!is.null(getPowerSim())){
    
    if(getPowerSim()$ESTmethod=="SP") mypars = c( "beta", "AUC", "TPR", "FPR", "PPV", "NPV")
    else mypars = c( "AUC", "TPR", "FPR", "PPV", "NPV")
    
    
    mypars = mypars[mypars!=substr(input$parameter, 1,3)]
    if(length(mypars)>0){
      printResultPlot(getPowerSim(), pars= mypars, useLogit = FALSE)
    }else{ return()}
    
  }
})

  
  
  
})
