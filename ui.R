
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
shinyUI(fluidPage(

  # Application title
  titlePanel("Comparison of GCMs models"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("year", "Years to include:",
                         c("2050" = "2050",
                           "2070" = "2070"), inline = TRUE, selected = c("2050")),
      checkboxGroupInput("rcp", "Select RCPs to inlcude:",
                         c("RCP 2.6" = "RCP 2.6",
                           "RCP 4.5" = "RCP 4.5",
                           "RCP 6.0" = "RCP 6.0",
                           "RCP 8.5" = "RCP 8.5"), inline = TRUE, selected = c("RCP 4.5")),
      checkboxGroupInput("all.models", "Select Models to inlcude:",
                         c("AC", "BC", "CC", "CE", "CN", "GF", "GD", "GS", "HD", "HG", "HE", "IN", "IP", "MI", "MR", "MC", "MP", "MG","NO"), inline = TRUE, selected = c("BC", "CC", "CE")),
      checkboxGroupInput("selected.bio", "Bioclim Variables to include:",
                         c("Bio1" = 1, "Bio2" = 2, "Bio3" = 3, "Bio4" = 4, "Bio5" = 5, "Bio6" =6, "Bio7" = 7, "Bio8" = 8, "Bio9" = 9, "Bio10" = 10, "Bio11" = 11, "Bio12" = 12, "Bio13" = 13, "Bio14" = 14, "Bio15" = 15, "Bio16" = 16, "Bio17" = 17, "Bio18" = 18,"Bio19" = 19), inline = TRUE, selected = c(1, 5)),
      #actionButton("do", "Download layers"),
      radioButtons("type", "Choose how will you define your extent:",
                   c("Numeric input" = "num",
                     "By clicking" = "click",
                     "By choosing counties" = "con")),
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
                      c("World", "Antigua and Barbuda", "Algeria", "Azerbaijan", "Albania", 
                        "Armenia", "Angola", "American Samoa", "Argentina", "Australia", 
                        "Bahrain", "Barbados", "Bermuda", "Bahamas", "Bangladesh", "Belize", 
                        "Bosnia and Herzegovina", "Bolivia", "Burma", "Benin", "Solomon Islands", 
                        "Brazil", "Bulgaria", "Brunei Darussalam", "Canada", "Cambodia", 
                        "Sri Lanka", "Congo", "Democratic Republic of the Congo", "Burundi", 
                        "China", "Afghanistan", "Bhutan", "Chile", "Cayman Islands", 
                        "Cameroon", "Chad", "Comoros", "Colombia", "Costa Rica", "Central African Republic", 
                        "Cuba", "Cape Verde", "Cook Islands", "Cyprus", "Denmark", "Djibouti", 
                        "Dominica", "Dominican Republic", "Ecuador", "Egypt", "Ireland", 
                        "Equatorial Guinea", "Estonia", "Eritrea", "El Salvador", "Ethiopia", 
                        "Austria", "Czech Republic", "French Guiana", "Finland", "Fiji", 
                        "Falkland Islands (Malvinas)", "Micronesia, Federated States of", 
                        "French Polynesia", "France", "Gambia", "Gabon", "Georgia", "Ghana", 
                        "Grenada", "Greenland", "Germany", "Guam", "Greece", "Guatemala", 
                        "Guinea", "Guyana", "Haiti", "Honduras", "Croatia", "Hungary", 
                        "Iceland", "India", "Iran (Islamic Republic of)", "Israel", "Italy", 
                        "Cote d'Ivoire", "Iraq", "Japan", "Jamaica", "Jordan", "Kenya", 
                        "Kyrgyzstan", "Korea, Democratic People's Republic of", "Kiribati", 
                        "Korea, Republic of", "Kuwait", "Kazakhstan", "Lao People's Democratic Republic", 
                        "Lebanon", "Latvia", "Belarus", "Lithuania", "Liberia", "Slovakia", 
                        "Liechtenstein", "Libyan Arab Jamahiriya", "Madagascar", "Martinique", 
                        "Mongolia", "Montserrat", "The former Yugoslav Republic of Macedonia", 
                        "Mali", "Morocco", "Mauritius", "Mauritania", "Malta", "Oman", 
                        "Maldives", "Mexico", "Malaysia", "Mozambique", "Malawi", "New Caledonia", 
                        "Niue", "Niger", "Aruba", "Anguilla", "Belgium", "Hong Kong", 
                        "Northern Mariana Islands", "Faroe Islands", "Andorra", "Gibraltar", 
                        "Isle of Man", "Luxembourg", "Macau", "Monaco", "Palestine", 
                        "Montenegro", "Mayotte", "Aaland Islands", "Norfolk Island", 
                        "Cocos (Keeling) Islands", "Antarctica", "Bouvet Island", "French Southern and Antarctic Lands", 
                        "Heard Island and McDonald Islands", "British Indian Ocean Territory", 
                        "Christmas Island", "United States Minor Outlying Islands", "Vanuatu", 
                        "Nigeria", "Netherlands", "Norway", "Nepal", "Nauru", "Suriname", 
                        "Nicaragua", "New Zealand", "Paraguay", "Peru", "Pakistan", "Poland", 
                        "Panama", "Portugal", "Papua New Guinea", "Guinea-Bissau", "Qatar", 
                        "Reunion", "Romania", "Republic of Moldova", "Philippines", "Puerto Rico", 
                        "Russia", "Rwanda", "Saudi Arabia", "Saint Kitts and Nevis", 
                        "Seychelles", "South Africa", "Lesotho", "Botswana", "Senegal", 
                        "Slovenia", "Sierra Leone", "Singapore", "Somalia", "Spain", 
                        "Saint Lucia", "Sudan", "Sweden", "Syrian Arab Republic", "Switzerland", 
                        "Trinidad and Tobago", "Thailand", "Tajikistan", "Tokelau", "Tonga", 
                        "Togo", "Sao Tome and Principe", "Tunisia", "Turkey", "Tuvalu", 
                        "Turkmenistan", "United Republic of Tanzania", "Uganda", "United Kingdom", 
                        "Ukraine", "United States", "Burkina Faso", "Uruguay", "Uzbekistan", 
                        "Saint Vincent and the Grenadines", "Venezuela", "British Virgin Islands", 
                        "Viet Nam", "United States Virgin Islands", "Namibia", "Wallis and Futuna Islands", 
                        "Samoa", "Swaziland", "Yemen", "Zambia", "Zimbabwe", "Indonesia", 
                        "Guadeloupe", "Netherlands Antilles", "United Arab Emirates", 
                        "Timor-Leste", "Pitcairn Islands", "Palau", "Marshall Islands", 
                        "Saint Pierre and Miquelon", "Saint Helena", "San Marino", "Turks and Caicos Islands", 
                        "Western Sahara", "Serbia", "Holy See (Vatican City)", "Svalbard", 
                        "Saint Martin", "Saint Barthelemy", "Guernsey", "Jersey", "South Georgia South Sandwich Islands", 
                        "Taiwan"), multiple = TRUE, selected = "World") 
                      #,submitButton("Update View", icon("refresh"))
                      )
      ),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot", click = "plot_click", brush = "plot_brush"),
      h2("You have selected"),
      verbatimTextOutput("visFun"),
      verbatimTextOutput("info"),
      dataTableOutput("table")
    )
  )
))
