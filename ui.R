

# Define UI for dataset viewer application
shinyUI(pageWithSidebarNavbar(
  includeHTML("www/navbar.html"),

  # Application title
 headerPanel(
    
    tags$body(
      h2("Power calculations for prognostic accuracy measures with survival data"), 
      h3(HTML("  &nbsp;&nbsp;<FONT COLOR=\"#f89406\"><em>Nested Case-Control study design</em></FONT>")) 
    ),
    windowTitle = "Power calculations for prognostic accuracy measures"
    
  ),
  
  sidebarPanel(
    includeHTML("www/navbar_setmethod.html"),
    wellPanel(
      selectInput("parameter", "Accuracy Measure", choices = c("AUC", "TPR(c)", "FPR(c)", "PPV(c)", "NPV(c)")),
      sliderInput("cutoff", "Cutpoint for marker (c)",min = -2, max = 2, value = 0,step = .05), #uiOutput("cutoffInput"), 
      p("note: cutoff value applies to all measures besides the AUC")
    ),
    wellPanel(
     selectInput("ESTmethod", "Estimation Method:", 
                 choices = c( "Non-parametric" = "NP", "Semi-parametric" = "SP")), 
     p("note: Non-parametric estimation methods will be added soon")
    ),
    wellPanel(p("Model Parameters:"),
              numericInput("S.0", label = HTML("Baseline Survival at time t<sub>0</sub>: S(t<sub>0</sub>|Y=0)"), value = round(exp(-.1), 3), min = .01, max = .99, step = .01),
              numericInput("t.0", label = HTML("time t<sub>0</sub>"), 1),
              numericInput("predict.time", "Future prediction time", 2)
    ),
    wellPanel(
              selectInput("censorType", "Type of Censoring", 
                                 choices = c("Percentage of observations" = "cens.perc", 
                                             "Maximum Time" = "maxtime")),
              uiOutput("censType")
             # sliderInput("cens.perc", label = "Percentage Censoring:", min = 0, max = 99, value = 20 )
      
      ),

    wellPanel( p("Hypothesis Test:"),

               htmlOutput("NullLine1"), 
               htmlOutput("NullLine2"), 
               uiOutput("nullInput"),

               numericInput("alpha", "Alpha level for test:", 0.05, min = 0.0, max = .99, step = .005),
               br(), br(),
               uiOutput("altInput")
    )
    #wellPanel(checkboxInput('useLogit', label="Use Logit Transform?", FALSE), 
    #          p(""))
    
  ),
  
  mainPanel( 
    tags$head(
      tags$script(src = 'https://c328740.ssl.cf1.rackcdn.com/mathjax/2.0-latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML', type = 'text/javascript') ,
        tags$title('Table Formatting'),
        tags$link(rel = 'stylesheet', type = 'text/css', href = 'stylesheet.css')
      ),
    tabsetPanel(
           tabPanel("Introduction", 
               includeRmd("introduction.Rmd"), br()
                    ),
           tabPanel("Instructions", 
               includeRmd("instructions.Rmd"), br()
                    ),
      tabPanel("Model Characteristics", 
               p("It may take a moment to display figures."),
               actionButton(inputId = "plotNow", "Generate Figures"),
               h3("Distribution of Survival and Y"),
               p(HTML("For survival time T, the figure on the left shows survival, S(t |Y) = Pr(T > t |Y), with respect to time for several different marker values. The vertical line highlights the prediction time. On the lower right, the marker distribution, along with cutpoint is shown. 
                 The upper right shows S(t|Y) by marker value when t equals the prediction time chosen.")),
               plotOutput("DistributionPlot",width="750px",height="350px"),
               div(class="row",
               #div(class="span2",  numericInput("lowertime", label = "Minimum time", value = 0 )),
               #div(class ="span1", br()),
               #div(class="span2",  numericInput("uppertime", label = "Maximum time", value = 20 ))),
                   div(class = "span1", br()),
                   div(class = "span4", numericInput("uppertime", label = "max time (for figure axis)", value = 20))),
               
               br(),
               htmlOutput("getEventRates"),
               htmlOutput("censoringPercentage"),
               br(),
               h3("Performance"),
               p(HTML("Given the inputs received, the table below shows the true value for each measure, including β, the coefficient in the Cox-Proportional hazards model. Also shown are the values of all measures under H<sub>0</sub>.
                 The figure below shows the true receiver operating characteristic (ROC) curve, along with the curve under the null hypothesis. In addition, curves for PPV(c) and NPV(c) are shown with respect to marker quantile. The chosen cutpoint is highlighted on each curve.")),
               div(class="row",
               div(class ="span1", br()),
               div(class="span4", tableOutput("TrueValuesTable")),
               div(class ="span1", br()),
               div(class="span6",  plotOutput("CurvesPlot", height = "500px")))
                             ),
           
     tabPanel("Set Sample Size", 
              wellPanel(
                numericInput("cohortN", "Cohort Sample Size", min = 25, max = 5000, value = 500),
                numericInput("nmatch", "Number of matched controls per case", min = 1, max = 10, value = 3)
              ),
              br(),
              htmlOutput("NCCsamplesizeALL"),
              br(),
              h5("By censoring status:"), 
              tableOutput("NCCsamplesize"), br(),
              #div(class = "row",
              #  div(class ="span1", br()),
              #  div(class = "span2",p("(Average) number of observations in:")),
              #  div(class = "span4", tableOutput("NCCsamplesize"))
              #    )
              HTML("<em>Note: To change the proportion of samples with observed/censored failure time in the cohort, change the censoring parameters to the left. </em>")
              ),
      
      tabPanel("Simulate Power", 

               p("This simulation takes many minutes to run (20+ minutes for cohort sample sizes > 500). Semiparametric estimation takes less computation time on average than nonparametric estimation. It is highly recommended to run this app from your local machine for larger simulations. See tab 'Instructions' for help doing this."), 
               wellPanel(
                 numericInput("mmSim", "# of Simulations to Run (enter a value between 10 and 1,000)", value = 100, step = 50, min = 10, max = 1000),
                 actionButton(inputId = "runSim", "Run Simulation")
               ),
               uiOutput("printEstMethod"), 
               uiOutput("printSampSize"), 
               uiOutput("printSampSizeNCC"), br(),
               uiOutput("PowerResultsHeader"), 
               plotOutput("simulationGraphTop",width="650px",height="180px"),
               br(), 
               h3("Other Measures:"),
               plotOutput("simulationGraphSub",width="650px",height="900px"),
                
               value = 2), id = 1
    ),
    
    ### show timer
    conditionalPanel("updateBusy() || $('html').hasClass('shiny-busy')",
                     id='progressIndicator',
                     "Calculation IN PROGRESS...",
                     div(id='progress',includeHTML("timer.js"))
    ),
    
    tags$head(tags$style(type="text/css",
                         '#progressIndicator {',
                         '  position: fixed; top: 110px; right: 8px; width: 200px; height: 50px;',
                         '  padding: 8px; border: 1px solid #CCC; border-radius: 8px;',
                         '}'
    ))
  )
))
