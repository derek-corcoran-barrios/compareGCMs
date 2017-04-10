
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(magrittr)
library(dplyr)
library(raster)
library(ccafs)
library(rasterVis)
library(maptools)
library(shiny)
library(rgeos)
library(sp)
library(rworldxtra)

shinyUI(fluidPage(
  
  # Application title
  titlePanel(div(img(src="sparc_logo.png", height=57,width=131), "Comparison of GCMs models")),
  titlePanel(h4("This App will assist you in exploring and comparing climate projections from 
                different Global Circulation Models (GCM) for projecting species distribution models")),
  titlePanel(h5("Make your choices in the left panel and find your results on the 'explore' and 'compare' tabs")),
  
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      h4("SET YOUR PARAMETERS"),
      h5("1. CHOOSE THE SCENARIOS AND MODELS TO EVALUATE:"),
      checkboxGroupInput("year", "- Years:",
                         c("2050" = "2050",
                           "2070" = "2070"), inline = TRUE, selected = c("2050")),
      checkboxGroupInput("rcp", "- RCPs:",
                         c("RCP 2.6" = "RCP 2.6",
                           "RCP 4.5" = "RCP 4.5",
                           "RCP 6.0" = "RCP 6.0",
                           "RCP 8.5" = "RCP 8.5"), inline = TRUE, selected = c("RCP 4.5")),
      checkboxGroupInput("all.models", "- Global Circulation Models (GCM): -find possible combinations in tab 'possible scenarios'-",
                         c("AC", "BC", "CC", "CE", "CN", "GF", "GD", "GS", "HD", "HG", "HE", "IN", "IP", "MI", "MR", "MC", "MP", "MG","NO"), inline = TRUE, selected = c("BC", "CC", "CE")),
      checkboxGroupInput("selected.bio", "- Bioclim Variables:",
                         c("Bio1" = 1, "Bio2" = 2, "Bio3" = 3, "Bio4" = 4, "Bio5" = 5, "Bio6" =6, "Bio7" = 7, "Bio8" = 8, "Bio9" = 9, "Bio10" = 10, "Bio11" = 11, "Bio12" = 12, "Bio13" = 13, "Bio14" = 14, "Bio15" = 15, "Bio16" = 16, "Bio17" = 17, "Bio18" = 18,"Bio19" = 19), inline = TRUE, selected = c(1, 5)),
      h5("\n", "2. ZOOM-IN TO YOUR AREA OF INTEREST:"),
      radioButtons("type", "(the map will refresh when you press 'Update extent')",
                   c("Select extent on the map" = "click",
                     "Choose counties" = "con",
                     "Write the bounding-box coordinates" = "num")),
      #submitButton("Update View", icon("refresh")),
      conditionalPanel(condition = "input.type == 'num'",
                       div(style="display:inline-block",numericInput("minlon", "Most Western longitude (negative for western hemisfere):", -180, min = -180, max = 180, step = 0.5)),
                       div(style="display:inline-block",numericInput("maxlon", "Most Eastern longitude (negative for western hemisfere):", 180, min = -180, max = 180, step = 0.5)),
                       numericInput("minlat", "Most Southern latitude (negative for southern hemisfere):", -90, min = -90, max = 90, step = 0.5),
                       numericInput("maxlat", "Most Northern latitude (negative for southern hemisfere):", 90, min = -90, max = 90, step = 0.5)
                       #,submitButton("Update View", icon("refresh"))
      ),
      conditionalPanel(condition = "input.type == 'con'",
                       selectInput("country", "Or select a counrty",
                                   c("World","Aruba", "Afghanistan", "Angola", "Anguilla", "Albania", "Aland", 
                                     "Andorra", "United Arab Emirates", "Argentina", "Armenia", "American Samoa", 
                                     "Antarctica", "Ashmore and Cartier Islands", "French Southern and Antarctic Lands", 
                                     "Antigua and Barbuda", "Australia", "Austria", "Azerbaijan", 
                                     "Burundi", "Belgium", "Benin", "Burkina Faso", "Bangladesh", 
                                     "Bulgaria", "Bahrain", "The Bahamas", "Bosnia and Herzegovina", 
                                     "Saint Barthelemy", "Belarus", "Belize", "Bermuda", "Bolivia", 
                                     "Brazil", "Barbados", "Brunei", "Bhutan", "Botswana", "Central African Republic", 
                                     "Canada", "Switzerland", "Chile", "China", "Ivory Coast", "Clipperton Island", 
                                     "Cameroon", "Cyprus No Mans Area", "Democratic Republic of the Congo", 
                                     "Republic of the Congo", "Cook Islands", "Colombia", "Comoros", 
                                     "Cape Verde", "Costa Rica", "Coral Sea Islands", "Cuba", "Curacao", 
                                     "Cayman Islands", "Northern Cyprus", "Cyprus", "Czech Republic", 
                                     "Germany", "Djibouti", "Dominica", "Denmark", "Dominican Republic", 
                                     "Algeria", "Ecuador", "Egypt", "Eritrea", "Dhekelia Sovereign Base Area", 
                                     "Spain", "Estonia", "Ethiopia", "Finland", "Fiji", "Falkland Islands", 
                                     "France", "Faroe Islands", "Federated States of Micronesia", 
                                     "Gabon", "Gaza", "United Kingdom", "Georgia", "Guernsey", "Ghana", 
                                     "Gibraltar", "Guinea", "Gambia", "Guinea Bissau", "Equatorial Guinea", 
                                     "Greece", "Grenada", "Greenland", "Guatemala", "Guam", "Guyana", 
                                     "Hong Kong S.A.R.", "Heard Island and McDonald Islands", "Honduras", 
                                     "Croatia", "Haiti", "Hungary", "Indonesia", "Isle of Man", "India", 
                                     "Indian Ocean Territories", "British Indian Ocean Territory", 
                                     "Ireland", "Iran", "Iraq", "Iceland", "Israel", "Italy", "Jamaica", 
                                     "Jersey", "Jordan", "Japan", "Baykonur Cosmodrome", "Siachen Glacier", 
                                     "Kazakhstan", "Kenya", "Kyrgyzstan", "Cambodia", "Kiribati", 
                                     "Saint Kitts and Nevis", "Korea No Mans Area", "South Korea", 
                                     "Kosovo", "Kuwait", "Laos", "Lebanon", "Liberia", "Libya", "Saint Lucia", 
                                     "Liechtenstein", "Sri Lanka", "Lesotho", "Lithuania", "Luxembourg", 
                                     "Latvia", "Macau S.A.R", "Saint Martin", "Morocco", "Monaco", 
                                     "Moldova", "Madagascar", "Maldives", "Mexico", "Marshall Islands", 
                                     "Macedonia", "Mali", "Malta", "Myanmar", "Montenegro", "Mongolia", 
                                     "Northern Mariana Islands", "Mozambique", "Mauritania", "Montserrat", 
                                     "Mauritius", "Malawi", "Malaysia", "Namibia", "New Caledonia", 
                                     "Niger", "Norfolk Island", "Nigeria", "Nicaragua", "Niue", "Netherlands", 
                                     "Norway", "Nepal", "Nauru", "New Zealand", "Oman", "Pakistan", 
                                     "Panama", "Pitcairn Islands", "Peru", "Philippines", "Palau", 
                                     "Papua New Guinea", "Poland", "Puerto Rico", "North Korea", "Portugal", 
                                     "Paraguay", "French Polynesia", "Qatar", "Romania", "Russia", 
                                     "Rwanda", "Western Sahara", "Saudi Arabia", "Sudan", "South Sudan", 
                                     "Senegal", "Singapore", "South Georgia and South Sandwich Islands", 
                                     "Saint Helena", "Solomon Islands", "Sierra Leone", "El Salvador", 
                                     "San Marino", "Somaliland", "Somalia", "Saint Pierre and Miquelon", 
                                     "Republic of Serbia", "Sao Tome and Principe", "Suriname", "Slovakia", 
                                     "Slovenia", "Sweden", "Swaziland", "Sint Maarten", "Seychelles", 
                                     "Syria", "Turks and Caicos Islands", "Chad", "Togo", "Thailand", 
                                     "Tajikistan", "Turkmenistan", "East Timor", "Tonga", "Trinidad and Tobago", 
                                     "Tunisia", "Turkey", "Tuvalu", "Taiwan", "United Republic of Tanzania", 
                                     "Uganda", "Ukraine", "United States Minor Outlying Islands", 
                                     "Uruguay", "United States of America", "US Naval Base Guantanamo Bay", 
                                     "Uzbekistan", "Vatican", "Saint Vincent and the Grenadines", 
                                     "Venezuela", "British Virgin Islands", "United States Virgin Islands", 
                                     "Vietnam", "Vanuatu", "West Bank", "Wallis and Futuna", "Akrotiri Sovereign Base Area", 
                                     "Samoa", "Yemen", "South Africa", "Zambia", "Zimbabwe"), multiple = TRUE, selected = "World") 
                       #,submitButton("Update View", icon("refresh"))
      ),
      h5("\n", "3. PRESS THE BUTTON TO PROCEED:"),
      actionButton("go", "Update variables and extent")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(tabsetPanel(
      tabPanel("Main",plotOutput("distPlot", click = "plot_click", brush = "plot_brush", height="800px"),
               h2("You have selected"),
               verbatimTextOutput("visFun"),
               verbatimTextOutput("info")),
               # dataTableOutput("table")),
      tabPanel("EXPLORE",
               h4("These are the layers and GCMs you have chosen. Here you can inspect them and decide
                  if any of them does not have enough quality to be used in the comparison"),
               uiOutput("explore.plot")),
               # plotOutput("explore.plot")),
      tabPanel("COMPARE",
               h3("Table summarizing comparisons with ensembles"),
               h5("This table contains the results of the comparison exercise. You might filter different
                  scenarios using the boxes below each column. The values in each row result from comparing
                  chosen GCMs in the same scenario, and are obtained after summarizing how much
                  deviate from the average prediction accumulatively all pixels"),
               h5("- Values close to 'cero' -for a variable or in total for a GCM at a given scenario- indicate that 
                  the values in that model are equal to the average prediction among all models compared"),
               h5("- -1 indicates that the values in that the model has the most different and negative
                  predictions (i.e.: colder temperatures, rarer precipitations)"),
               h5("- +1 indicates that the values in that the model has the most different and positive
                  predictions (i.e.: warmer temperatures, more precipitations)"),
               dataTableOutput("table"),
               h3("Plots of comparisons with ensembles"),
               h5("Red colors/positive values indicate positive deviations from average prediction
                  among compared GCMs (i.e.: more temperature/precipitation), blue colors(negative
                  values indicate negative deviations"),
               uiOutput("compare.plot"), h5("If you have detected patterns that you do not like in any GCM, 
                  you might consider repeating the analysis after removing those GCM,
                                            and working with a comparison which is not 'contaminated' by 'bad' models")),
      tabPanel("Possible scenarios",
               dataTableOutput("Possible"))
               
               # plotOutput("compare.plot"))

      # tabPanel("plots", uiOutput('myBiotabs'))
    ))
  )
))

