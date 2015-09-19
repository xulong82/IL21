library(shiny)

shinyUI(
  navbarPage("Lupus Project",
    tabPanel("Expression",
    fluidPage(
    p(strong("Figure 1:"), "Gene expression profiles (TPM) in Naive, Act, and Il21-Act cells."),
    hr(),
    sidebarLayout(
      sidebarPanel(
        textInput("gene", label = h3("Input gene:"), value = "Bcl6")
      ),
      mainPanel(
        plotOutput("plot"),
	      tableOutput("summary")
      )
    )
  )),
	tabPanel("Pathway",
	  fluidRow(
	    column(2,
	      radioButtons("cell", "Choose cell", choices = list("N" = 1, "ACT" = 2, "ACT IL21" = 3), selected = 3),
        radioButtons("layout", "Choose graph layout", choices = list("Sphere" = 1, "Circle" = 2, "Community" = 3), selected = 1)
	    ),
	    column(10,
	      plotOutput("igraph")
	    )
	)),
	tabPanel("Documents",
	  h3("Useful documents"),
	  hr(),
	  p("Complete gene expression in TPM unit"),
	  downloadButton('downloadData', 'Download'),
	  hr(),
	  p("Manuscript: Supplementary Table 1"),
	  downloadButton('downloadDataM1', 'Download'),
	  hr(),
	  p("Manuscript: Supplementary Table 2"),
	  downloadButton('downloadDataM2', 'Download'),
	  hr(),
	  h5("Methods on the computational side"),
	  p(a(href="profile.html", "Identify signature genes")),
	  p(a(href="immgen.html", "Integrate with IMMGEN dataset")),
	  p(a(href="liu.html", "Integrate with Liu et al.")),
	  hr(),
    p(a("GitHub Repository", href = "https://github.com/xulong82/Lupus"))
	)
))