PowerSAM_ncc
==============

Simulate power for prognostic biomarker validation studies with time-to-event data under nested case control design. 
The application can be found [here](http://mdbrown.shinyapps.io/PowerSAM_ncc/), or, if you would rather run it locally on your machine use:

```r
#install package shiny
if (!require("shiny")) install.packages("shiny")

#install devtools so we can install shinyIncubator
if(!require("devtools")) install.packages("devtools")
if(!require("shinyIncubator")) devtools::install_github("shiny-incubator", "rstudio")


#now run the app locally using
runGitHub("PowerSAM_ncc", "mdbrown")


```