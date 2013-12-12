includeRmd <- function(path){
  if (!require(knitr))
    stop("knitr package is not installed")
  if (!require(markdown))
    stop("Markdown package is not installed")
  shiny:::dependsOnFile(path)
  contents = paste(readLines(path, warn = FALSE), collapse = '\n')
  html <- knitr::knit2html(text = contents, fragment.only = TRUE)
  Encoding(html) <- 'UTF-8'
  return(HTML(html))
}

myheaderPanel <- function(title, windowTitle=title) {    
  tagList(
    tags$head(tags$title(windowTitle)),
    div(class="span12", style="padding: 10px 0px;",
        h3(title)
    )
  )
}



pageWithSidebarNavbar <- function(navbar, headerPanel, sidebarPanel, mainPanel) {
  
  bootstrapPage(
    # basic application container divs
    navbar,
    div(
      class="container-fluid", 
      div(class="row-fluid", 
          headerPanel
      ),
      div(class="row-fluid", 
          sidebarPanel, 
          mainPanel
      )
    )
  )
}