library(shiny)

shinyUI(
  fluidPage(
    # CSS theme
    #theme = "css/bootstrap.css",
    
    # Application title
    titlePanel("Mathematical Models of Infectious Diseases in Human and Animals"),
    
    # Sidebar
    sidebarLayout(
      sidebarPanel(width = 4,
        h4("Models in SIR:"),
        
        selectInput(inputId = "modtype",
                    label = "",
                    choices = c("SIR", "SIR + 2 Age Class", "SIR + 2 Type Imports",
                    "SIR + Constant Additive","SIR + Birth Death","SIR + Carrier State",
                    "SIR + DemogStoch","SIR + Induced Mortality","SIR + Partial Immunity",
                    "SIR + Sinusoidal Births","SIR + Sinusoidal Forcing","SIR + Tau Leap",
                    "SIR + Vector"),selected="SIR"),
        br(),
        
        h4("Parameters:"),
        
 ######################## SIR   #####################################
        
      conditionalPanel(
          condition = "input.modtype == 'SIR'",
          numericInput(inputId = "beta",label = "Transmission rate:",
                             min = 0,max = 10,value = 0.5),br(),  ####### for beta
          numericInput(inputId = "gamma",label = "Removal rate:",
                             min = 0,max = 10,value = 0.5),br(),br(), ####### for gamma
        
         h4("Initial Proportion:"),
        
         numericInput(inputId = "S",label = "Proportion of Suspcetible:",
                      min = 0, max = 1,value = 0.01,step=0.1),     
         numericInput(inputId = "I",label = "Proportion of Infectious:",
                     min = 0, max = 1,value = 0.01, step=0.1),
         numericInput(inputId = "R",label = "Proportion of Recovery:",
                      min = 0, max = 1,value = 0.01, step=0.1),br(),    
         h4("Other parameters:"),
         numericInput(inputId="time",label="Time Sequence:",
                      min=1,max=1e10,value=365,step=365)),
 
    ##################### SIR + 2 Age Class  #####################
        
    conditionalPanel(
          condition = "input.modtype == 'SIR + 2 Age Class'",
          numericInput(inputId = "betaCC",label = "Transmission rate for Children:",
                       min = 0,max = 1000,value = 100),br(),
          numericInput(inputId = "betaAA",label = "Transmission rate for Adult:",
                       min = 0,max = 1000,value = 20),br(),
          numericInput(inputId = "betaCA",label = "Transmission rate for Children to Adult:",
                       min = 0,max = 1000,value = 10),br(),
          numericInput(inputId = "betaAC",label = "Transmission rate for Adult to Children:",
                       min = 0,max = 1000,value = 10),br(),
          numericInput(inputId = "gamma",label = "Removal rate:",
                       min = 0,max = 100,value = 10),br(),
          numericInput(inputId = "lC",label = "Children muture rate:",
                       min = 0,max = 1,value = 0.067),br(),
          numericInput(inputId = "muC",label = "Children death rate:",
                       min = 0,max = 1,value = 0.01),br(),
          numericInput(inputId = "muA",label = "Adult death rate:",
                       min = 0,max = 1,value = 0.0167),br(),br(),
          
          h4("Initial Proportion:"),
          
          numericInput(inputId = "SC",
                       label = "Proportion of Suspcetible in Children:",
                       min = 0, max = 1,value = 0.1,step=0.1),     
          numericInput(inputId = "IC",
                       label = "Proportion of Infectious in Children:",
                       min = 0, max = 1,value = 0.001, step=0.1),
          numericInput(inputId = "SA",
                       label = "Proportion of Infectious in Adult:",
                       min = 0, max = 1,value = 0.1, step=0.1),
          numericInput(inputId = "IA",
                       label = "Proportion of Infectious in Adult:",
                       min = 0, max = 1,value = 0.001, step=0.1),br(), 
   
          h4("Other parameters:"),
          
          numericInput(inputId="time",label="Time Sequence:",
                       min=1,max=1e10,value=365)),
 
    ##################### SIR + 2 Type Imports  ########################
    
    conditionalPanel(
      condition = "input.modtype == 'SIR + 2 Type Imports'",
      numericInput(inputId = "beta",label = "Transmission rate:",
                   min = 0,max = 10,value = 0.5),br(),  ####### for beta
      numericInput(inputId = "gamma",label = "Removal rate:",
                   min = 0,max = 10,value = 0.1),br(), ####### for gamma
      numericInput(inputId = "mu",label = "Per capital death rate:",
                   min = 0,max = 1,value = 5e-4),br(), ####### for mu
      numericInput(inputId = "epsilon",label = "Import via immigration rate:",
                   min = 0,max = 1,value = 2e-5),br(), ####### for epsilon
      numericInput(inputId = "delta",label = "Import via extenrnal infection rate:",
                   min = 0,max = 1,value = 0.01),br(),br(),####### for delta
    
      h4("Initial number:"),
      
      numericInput(inputId = "X",label = "Number of Suspcetible:",
                   min = 0, max = 1000,value = 5,step=100),     
      numericInput(inputId = "Y",label = "Number of Infectious:",
                   min = 0, max = 1000,value = 21, step=10),
      numericInput(inputId = "N",label = "Number of Total population:",
                   min = 0, max = 1e8,value = 200, step=100),br(),
           h4("Other parameters:"),
           numericInput(inputId="time",label="Time Sequence:",
                        min=1,max=1e10,value=365,step=365)),
 
    #####################   SIR + Constant Additive   ####################
    
    conditionalPanel(
      condition = "input.modtype == 'SIR + Constant Additive'",
      numericInput(inputId = "beta",label = "Transmission rate:",
                   min = 0,max = 10,value = 0.5),br(),  ####### for beta
      numericInput(inputId = "gamma",label = "Removal rate:",
                   min = 0,max = 10,value = 0.1),br(), ####### for gamma
      numericInput(inputId = "mu",label = "Per capital death rate:",
                   min = 0,max = 1,value = 1/(50*365)),br(), ####### for mu
      numericInput(inputId = "noise",label = " Noise in transmission rate:",
                   min = 0,max = 100,value = 10),br(), ####### for noise
      numericInput(inputId = "N",label = "Constant population size:",
                   min = 0,max = 10,value = 0.5),br(), ####### for num
      
      h4("Initial number:"),
      
      numericInput(inputId = "X",label = "Number of Suspcetible:",
                   min = 0, max = 1e10,value = 1e5,step=100),     
      numericInput(inputId = "Y",label = "Number of Infectious:",
                   min = 0, max = 1e10,value = 500, step=100),br(),
      
      h4("Other parameters:"),
      numericInput(inputId="time",label="Time Sequence:",
                   min=1,max=1e10,value=365,step=365)),
    
   ######################## SIR + Birth Death ####################
   
   conditionalPanel(
     condition = "input.modtype == 'SIR + Birth Death'",
     numericInput(inputId = "beta",label = "Transmission rate:",
                  min = 0,max = 10,value = 1.35),br(),  ####### for beta
     numericInput(inputId = "gamma",label = "Removal rate:",
                  min = 0,max = 10,value = 1/7),br(), ####### for gamma
     numericInput(inputId = "mu",label = "Per capital death rate:",
                  min = 0,max = 1,value = 1/(70*365)),br(),br(), ####### for mu
     
     h4("Initial proportion:"),
     
     numericInput(inputId = "S",label = "Proportion of Suspcetible:",
                  min = 0, max = 1,value = 0.1,step=0.1),     
     numericInput(inputId = "I",label = "Proportion of Infectious:",
                  min = 0, max = 1,value = 1e-3, step=0.1),br(),
       h4("Other parameters:"),
       numericInput(inputId="time",label="Time Sequence:",
                    min=1,max=1e10,value=365,step=365)),
 
   ######################### SIR + Carrier State #########################
   
   conditionalPanel(
     condition = "input.modtype == 'SIR + Carrier State'",
     numericInput(inputId = "beta",label = "Transmission rate:",
                  min = 0,max = 10,value = 0.2),br(),  ####### for beta
     numericInput(inputId = "gamma",label = "Infectious-recovery rate:",
                  min = 0,max = 10,value = 0.1),br(), ####### for gamma
     numericInput(inputId = "mu",label = "Per capital death rate:",
                  min = 0,max = 10,value = 1/(50*365)),br(), ####### for mu
     numericInput(inputId = "Gamma",label = "Carrier-recovery rate:",
                  min = 0,max = 10,value = 0.001),br(), ####### for Gamma
     numericInput(inputId = "epsilon",label = "Proportion of reduction in transmission from carriers:",
                  min = 0,max = 10,value = 0.1),br(), ####### for mu
     numericInput(inputId = "rho",label = "Proportion of infectious to carriers:",
                  min = 0,max = 10,value = 0.4),br(),br(), ####### for mu
     
     h4("Initial proportion:"),
     
     numericInput(inputId = "S",label = "Proportion of Suspcetible:",
                  min = 0, max = 1,value = 0.1,step=0.1),     
     numericInput(inputId = "I",label = "Proportion of Infectious:",
                  min = 0, max = 1,value = 1e-4, step=0.1),
     numericInput(inputId = "C",label = "Proportion of Carriers:",
                  min = 0, max = 1,value = 1e-3, step=0.1),br(),
     h4("Other parameters:"),
     numericInput(inputId="time",label="Time Sequence:",
                  min=1,max=1e10,value=365,step=365)),
   
###########################   SIR + DemogStoch    ##########################

  conditionalPanel(
    condition = "input.modtype == 'SIR + DemogStoch'",
    numericInput(inputId = "beta",label = "Transmission rate:",
                 min = 0,max = 10,value = 1),br(),  ####### for beta
    numericInput(inputId = "gamma",label = "Infectious-recovery rate:",
                 min = 0,max = 10,value = 0.1),br(), ####### for gamma
    numericInput(inputId = "mu",label = "Per capital death rate:",
                 min = 0,max = 1,value = 5e-4),br(),br(), ####### for mu
    
    h4("Initial number:"),
    
    numericInput(inputId = "X",label = "Number of Suspcetible:",
                 min = 0, max = 1e10,value =500,step=100),     
    numericInput(inputId = "Y",label = "Number of Infectious:",
                 min = 0, max = 1e10,value = 25, step=100),
    numericInput(inputId = "N",label = "Number of total:",
                 min = 0, max = 1e10,value = 5e3, step=100),br(),

          h4("Other parameters:"),
          numericInput(inputId="time",label="Time Sequence:",
                       min=1,max=1e10,value=365,step=365)),

###################### SIR + Induced Mortality   #######################
  
  conditionalPanel(
    condition = "input.modtype == 'SIR + Induced Mortality'",
    numericInput(inputId = "beta",label = "Transmission rate:",
                 min = 0,max = 10,value = 0.5),br(),  ####### for beta
    numericInput(inputId = "gamma",label = "Infectious-recovery rate:",
                 min = 0,max = 10,value = 0.5),br(), ####### for gamma
    numericInput(inputId = "mu",label = "Per capital death rate:",
                 min = 0,max = 10,value = 0.5),br(), ####### for mu
    numericInput(inputId = "rho",label = "Probability of infected individual \n dies before recovery:",
                 min = 0,max = 1,value = 0.01,step=0.01),br(), ####### for rho
    numericInput(inputId = "nu",label = "Population level birth rate:",
                 min = 0,max = 1,value = 0.01,step=0.01),br(),br(), ####### for nu
    
    h4("Initial number:"),
    
    numericInput(inputId = "X",label = "Number of Suspcetible:",
                 min = 0, max = 1,value = 0.2,step=0.1),     
    numericInput(inputId = "Y",label = "Number of Infectious:",
                 min = 0, max = 1,value = 1e-4, step=0.01),
    numericInput(inputId = "Z",label = "Number of recovery:",
                 min = 0, max = 1,value = 0, step=0.01), br(),
    h4("Other parameters:"),
    numericInput(inputId="time",label="Time Sequence:",
                 min=1,max=1e9,value=365,step=365)
  ),
 
##################  SIR + Partial Immunity  ##################

conditionalPanel(
  condition = "input.modtype == 'SIR + Partial Immunity'",
  numericInput(inputId = "beta1",label = "Transmission rate for strain 1:",
               min = 0,max = 10,value = 260/365),br(),  ####### for beta1
  numericInput(inputId = "beta2",label = "Transmission rate for strain 2:",
               min = 0,max = 10,value = 520/365),br(),  ####### for beta2
  numericInput(inputId = "gamma1",label = "Infectious-recovery rate for strain 1:",
               min = 0,max = 10,value = 1/7),br(), ####### for gamma1
  numericInput(inputId = "gamma2",label = "Infectious-recovery rate for strain 2:",
               min = 0,max = 10,value = 2/7),br(), ####### for gamma2
  numericInput(inputId = "mu",label = "Death rate:",
               min = 0,max = 1,value = 1/(70*365)),br(), ####### for mu
  numericInput(inputId = "v",label = "Birth rate:",
               min = 0,max = 10,value = 0.5),br(), ####### for v
  numericInput(inputId = "alpha1",label = "Modified susceptiblity for strain 1:",
               min = 0,max = 10,value = 0.5),br(), ####### for Gamma
  numericInput(inputId = "alpha2",label = "Modified susceptiblity for strain 2:",
               min = 0,max = 10,value = 0.4),br(), ####### for Gamma
  numericInput(inputId = "a1",label = "Modified transmission for strain 1:",
               min = 0,max = 10,value = 0.4),br(), ####### for a1
  numericInput(inputId = "a2",label = "Modified transmission for strain 2:",
               min = 0,max = 10,value = 0.5),br(),br(), ####### for a2
  
  h4("Initial proportion:"),
  
  numericInput(inputId = "NSS",label = "Proportion of NSS:",
               min = 0, max = 1,value = 0.1,step=0.1),     
  numericInput(inputId = "NRR",label = "Proportion of NRR:",
               min = 0, max = 1,value = 0.3798, step=0.1),
  numericInput(inputId = "NIS",label = "Proportion of NIS:",
               min = 0, max = 1,value = 1e-4, step=0.1),
  numericInput(inputId = "NRS",label = "Proportion of NRS:",
               min = 0, max = 1,value = 0.02, step=0.1),
  numericInput(inputId = "NRI",label = "Proportion of NRI:",
               min = 0, max = 1,value = 0, step=0.1),
  numericInput(inputId = "NSI",label = "Proportion of NSI:",
               min = 0, max = 1,value = 1e-4, step=0.1),
  numericInput(inputId = "NSR",label = "Proportion of NSR:",
               min = 0, max = 1,value = 0.5, step=0.1),
  numericInput(inputId = "NIR",label = "Proportion of NIR:",
               min = 0, max = 1,value = 0, step=0.1),br(),
  
  h4("Other parameters:"),
  numericInput(inputId="time",label="Time Sequence:",
               min=1,max=1e9,value=365,step=365)
),

########################## SIR + Sinusoidal Births ################

conditionalPanel(
  condition = "input.modtype == 'SIR + Sinusoidal Births'",
  numericInput(inputId = "beta",label = "Transmission rate:",
               min = 0,max = 10,value = 17/13,step=0.5),br(),  ####### for beta
  numericInput(inputId = "gamma",label = "Infectious-recovery rate:",
               min = 0,max = 1,value = 1/13,step=0.1),br(), ####### for gamma
  numericInput(inputId = "mu",label = "Per capital death rate:",
               min = 0,max = 1,value = 1/(50*365),step=0.01),br(), ####### for mu
  numericInput(inputId = "alpha0",label = "Mean birth rate:",
               min = 0,max = 1,value = 1/(50*365),step=0.01),br(), ####### for alph0
  numericInput(inputId = "alpha1",label = "Amplitude of sinuoidal forcing:",
               min = 0,max = 1,value = 0.25,step=0.01),br(), ####### for alph1
  numericInput(inputId = "w",label = "Frequency of the oscillations:",
               min = 0,max = 1,value = 2*pi/365,step=0.01),br(),br(), ####### for alph1
  
  h4("Initial proportion:"),
  
  numericInput(inputId = "S",label = "proportion of Suspcetible:",
               min = 0, max = 1,value = 1/17,step=0.01),     
  numericInput(inputId = "I",label = "Proportion of Infectious:",
               min = 0, max = 1,value = 1e-4, step=0.01)
),

######################### SIR + Sinusoidal Forcing ########################
  conditionalPanel(
    condition = "input.modtype == 'SIR + Sinusoidal Forcing'",
    numericInput(inputId = "beta0",label = "Mean transmission rate:",
                 min = 0,max = 10,value = 0.1,step=0.1),br(),               ####### for beta0
    numericInput(inputId = "beta0",label = "Mean transmission rate:",
                 min = 0,max = 10,value = 17/13,step=0.1),br(),               ####### for beta1
    numericInput(inputId = "gamma",label = "Infectious-recovery rate:",
                 min = 0,max = 1,value = 1/13,step=0.01),br(),               ####### for gamma
    numericInput(inputId = "mu",label = "Per capital death rate:",
                 min = 0,max = 1,value = 1/(50*365)),br(),                ####### for mu
    numericInput(inputId = "omega",label = "Frequency of the oscillations:",
                 min = 0,max = 1,value =2*pi/365,step=0.01),br(),br(), ####### for alph1
    
    h4("Initial proportion:"),
    
    numericInput(inputId = "S",label = "Proportion r of Suspcetible:",
                 min = 0, max = 1,value = 1/17,step=0.01),     
    numericInput(inputId = "I",label = "Proportion of Infectious:",
                 min = 0, max = 1,value = 1e-4, step=0.01)
  ),

########################## SIR + Tau Leap ################################

  conditionalPanel(
    condition = "input.modtype == 'SIR + Tau Leap'",
    numericInput(inputId = "beta",label = "Transmission rate:",
                 min = 0,max = 10,value = 1,step=0.1),br(),               ####### for beta
    numericInput(inputId = "gamma",label = "Infectious-recovery rate:",
                 min = 0,max = 10,value = 0.1,step=0.1),br(),               ####### for gamma
    numericInput(inputId = "mu",label = "Per capital death rate:",
                 min = 0,max = 10,value = 5e-4,step=0.1),br(),               ####### for mu
    numericInput(inputId = "tau",label = "Time step:",
                 min = 0,max = 100,value = 1,step=2),br(),       ####### for tau
    numericInput(inputId = "N",label = "Constant population size:",
                 min = 0,max = 1e10,value = 50,step=100),br(),br(), ####### for N
    
    h4("Initial number:"),
    
    numericInput(inputId = "X",label = "Number of Suspcetible:",
                 min = 0, max = 1e10,value = 5,step=100),     
    numericInput(inputId = "Y",label = "Number of Infectious:",
                 min = 0, max = 1e10,value = 1, step=100),  
    numericInput(inputId = "Z",label = "Number of Recovery:",
                 min = 0, max = 1e10,value = 44, step=100),br(),
    h4("Other parameters:"),
    numericInput(inputId="time",label="Time Sequence:",
                 min=1,max=1e10,value=365,step=365)
    ),


############################# SIR + Vector ##############################
  
  conditionalPanel(
    condition = "input.modtype == 'SIR + Vector'",
    numericInput(inputId = "muH",label = "Human mortality rate:",
                 min = 0,max = 1,value = 5.5e-5,step=0.001),br(),          ####### for muH
    numericInput(inputId = "muM",label = "Mosquito mortality rate:",
                 min = 0,max = 10,value = 0.143,step=0.01),br(),          ####### for muM
    numericInput(inputId = "vH",label = "Human birth rate:",
                 min = 0,max = 1,value = 5.5e-2,step=0.01),br(),          ####### for VH
    numericInput(inputId = "vM",label = "Mosquito birth rate:",
                 min = 0,max = 10e7,value = 1.44e3),br(),          ####### for VM
    numericInput(inputId = "betaHM",label = "Transmission Pobability of \n humans-mosquitos:",
                 min = 0,max = 10,value = 0.5,step=0.1),br(),    ####### for betaHM
    numericInput(inputId = "betaMH",label = "Transmission Pobability of \n mosquitos-humans:",
                 min = 0,max = 10,value = 0.8,step=0.1),br(),    ####### for betaMH
    numericInput(inputId = "gamma",label = "Infectious-recovery rate:",
                 min = 0,max = 10,value = 0.033,step=0.01),br(), ####### for gamma
    numericInput(inputId = "r",label = "Bitten rate:",
                 min = 0,max = 10,value = 0.5/1e3),br(),br(), ####### for mu
    
    h4("Initial number:"),
    
    numericInput(inputId = "XH",label = "Number of Suspcetible for Human:",
                 min = 0, max = 1e10,value = 1e3,step=1e2), 
    numericInput(inputId = "XM",label = "Number of Suspcetible for Mosquitos:",
                 min = 0, max = 1e10,value = 1e4,step=100), 
    numericInput(inputId = "YH",label = "Proportion of Infectious for Human:",
                 min = 0, max = 100,value = 1, step=1),
    numericInput(inputId = "YM",label = "Proportion of Infectious for Mosquitos:",
                 min = 0, max = 100,value = 1, step=2),br(),
    h4("Other parameters:"),
    numericInput(inputId="time",label="Time Sequence:",
                 min=1,max=1e10,value=365,step=365)
  )),
      
      mainPanel(
        tabsetPanel(
          tabPanel(
            title = "Plot",
            
            wellPanel(plotOutput("plot")),
            
            wellPanel(uiOutput("timeline"))
          ),     
          
          tabPanel("Data",
                   h4("Model Data"),
                   dataTableOutput("table"),
                   br(),br(),
                   downloadButton(outputId = "dlData",
                                  label = "Download Data")
          ),     
          
          # About panel
          tabPanel(
            title = "About",
            
            tags$hr(),
            
            p(strong("Author:"), " Spatial-R"),
            
            p(strong("Email:"), a("Bing Zhang", 
                                    href = "zhangbing4502431@outlook.com",
                                    target = "_blank")),
            
            p(strong("Website:"), a("http://spatial-r.github.io/", 
                                    href = "http://spatial-r.github.io/",
                                    target = "_blank")),
            
            p(strong("Source code:"), 
              a("GitHub",
                href = "https://github.com/Spatial-R/Dynamic-Model-for-Infectious-Disease",
                target = "_blank")),
            
            p(strong("Created with:"), 
              a("RStudio,",
                href = "http://www.rstudio.com/",
                target = "_blank"),
              a("Shiny",
                href = "http://shiny.rstudio.com",
                target = "_blank"),
              " and ",
              a("EpiDynamic.",
                href="https://github.com/oswaldosantos/EpiDynamics",
              target = "_blank")),
            
            p(strong("License:"), 
              a("GPL v3",
                href = "http://www.gnu.org/copyleft/gpl.html",
                target = "_blank")),
            
            tags$hr(),
            
            h4("Recommended links:"),
            
            p(HTML('<ol>'),
              HTML('<li>'), a("Physicists model how we form opinions",
                              href = "http://phys.org/news127385810.html",
                              target = "_blank"), HTML('</li>'),
              HTML('<li>'), a("Minority Rules: Scientists Discover Tipping Point for the Spread of Ideas",
                              href = "http://news.rpi.edu/luwakkey/2902",
                              target = "_blank"), HTML('</li>'),
              HTML('<li>'), p(a("Why Internet Rumors Spread So Quickly",
                                href = "http://www.youtube.com/watch?v=g8GKJ1GwFvg",
                                target = "_blank"), 
                              "(link to video below)"), HTML('</li>'),
              HTML('</ol>'))
          )
        )
      )
    )
  )
)

