
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

if (!require("pacman")) install.packages("pacman")
pacman::p_load(magrittr, dplyr, raster, ccafs, rasterVis, maptools, shiny, rgeos, sp, viridis, ggplot2, ggforce, rworldxtra, matrixStats)

shinyUI(fluidPage(
  
  # Application title
  titlePanel(div(img(src="sparc_logo.png", height=100,width=250), "GCMcompareR: Comparison of GCMs models")),
  # titlePanel("Comparison of GCMs models"),
  titlePanel(h4("This App will assist you in exploring and comparing climate projections from 
                different Global Circulation Models (GCM) for projecting species distribution models")),
  titlePanel(h5("Make your choices in the left panel and find your results on the 'explore' and 'compare' tabs")),
  
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      h4("SET YOUR PARAMETERS"),
      h5("1. SCENARIOS:"),
      radioButtons("year", "- Years:",
                         c("2050" = "2050",
                           "2070" = "2070"), inline = TRUE, selected = c("2070")
                         ),
      radioButtons("rcp", "- RCPs:",
                         c("RCP 2.6" = "26",
                           "RCP 4.5" = "45",
                           "RCP 6.0" = "60",
                           "RCP 8.5" = "85"), inline = TRUE, selected = c("45")
                         ),
      checkboxGroupInput("gcm", "- Global Circulation Models (GCM): -find possible combinations in tab 'possible scenarios'-",
                         c("access1_0","bcc_csm1_1","bnu_esm","cccma_canesm2","cesm1_cam5_1_fv2","cnrm_cm5","csiro_mk3_6_0","fio_esm","gfdl_cm3","gfdl_esm2g","gfdl_esm2m","giss_e2_h","giss_e2_r","HadGEM2_AO","inm_cm4","ipsl_cm5a_lr","ipsl_cm5a_mr","lasg_fgoals_g2","miroc_esm","miroc_esm_chem","miroc_miroc5","mohc_hadgem2_cc","mohc_hadgem2_es","mpi_esm_lr","mpi_esm_mr","mri_cgcm3","ncar_ccsm4","ncc_noresm1_m"), 
                         inline = TRUE, 
                         selected = c("access1_0","bcc_csm1_1","mpi_esm_lr","gfdl_esm2g","giss_e2_h","mpi_esm_mr","ncar_ccsm4","ncc_noresm1_m","mohc_hadgem2_cc","cnrm_cm5","gfdl_cm3","mohc_hadgem2_es","lasg_fgoals_g2","cesm1_cam5_1_fv2")
                         ),
      radioButtons("analysistype", "Choose one of these analysis:",
                   c("Temperature vs Precipitation comparison" = "tVSp",
                     "Custom comparison" = "custom")
                   ),
      conditionalPanel(condition = "input.analysistype == 'custom'",
                       checkboxGroupInput("selected.bio", "",
                                         c("Bio1" = 1, "Bio2" = 2, "Bio3" = 3, "Bio4" = 4, "Bio5" = 5, "Bio6" =6, "Bio7" = 7, "Bio8" = 8, "Bio9" = 9, "Bio10" = 10, "Bio11" = 11, "Bio12" = 12, "Bio13" = 13, "Bio14" = 14, "Bio15" = 15, "Bio16" = 16, "Bio17" = 17, "Bio18" = 18,"Bio19" = 19),
                                         inline = TRUE,
                                         selected = c(1,5,6,8,9,10,11, 12,13,14,16,17,18,19))
                       ),
                          
      h5("\n", "2. STUDY AREA"),
      radioButtons("type", "(the map will refresh when you press 'Update extent')",
                   c("Select extent on the map" = "click",
                     "Choose counties" = "con",
                     "Write the bounding-box coordinates" = "num")),
        conditionalPanel(condition = "input.type == 'num'",
                       div(style="display:inline-block",numericInput("minlon", "Most Western longitude (negative for western hemisfere):", -180, min = -180, max = 180, step = 0.5)),
                       div(style="display:inline-block",numericInput("maxlon", "Most Eastern longitude (negative for western hemisfere):", 180, min = -180, max = 180, step = 0.5)),
                       numericInput("minlat", "Most Southern latitude (negative for southern hemisfere):", -90, min = -90, max = 90, step = 0.5),
                       numericInput("maxlat", "Most Northern latitude (negative for southern hemisfere):", 90, min = -90, max = 90, step = 0.5)
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
                       # ,submitButton("Update View", icon("refresh"))
                       ),
      actionButton("updateExtent", "Update Extent", icon("refresh")),
      
      h5("\n", "3. COMPARE GCMS!"),
      actionButton("go", "Compare")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(type="pills",
            tabPanel("Main"
                     ,plotOutput("distPlot", click = "plot_click", brush = "plot_brush", height="800px")
                     # ,h2("You have selected")
                     # ,verbatimTextOutput("visFun")
                     # ,verbatimTextOutput("info")
                     #,dataTableOutput("table")
                  ),
            tabPanel("Possible scenarios",
                     dataTableOutput("Possible")),
            tabPanel("EXPLORE variables' plots",
                     h4("These are the layers and GCMs you have chosen as downloaded."),
                     h5("(this takes a minute to print, especially if you selected a large area)"),
                     uiOutput("GCMPlots")),

            tabPanel("COMPARE - Results",
                     h3("Comparisons result"),
                     h5("The table is sorted according to the simmilarity of each model to the ensemble predictions (lower dist,
                        which is, the distance from the center (coordinates 0,0) to the position of the model)."),
                     h5("It can be sorted according to temperature or precipitation, as desired"),
                     dataTableOutput("table"),
                     h3("Figure for Precipitation and temperature comparisons"),
                     h5("This figure represents together temperature and precipitation comparisons. Models closer to the center
                        are those more similar to the ensemble at both characteristics."),
                     plotOutput("pVSt.plot", width=600, height=600),
                     ## GCM CHOICE
                     h4("CHOOSING GCMS: "),
                     radioButtons("GCMchoice", "Pick a method to select GCMs:",
                                  c("Best bracketing" = "dif",
                                  "Most similar to ensemble" = "sim"),
                                  inline = T, 
                                  selected = "sim"),
                     conditionalPanel(condition = "input.GCMchoice == 'dif'",
                                      selectInput("hot_dry", "1. Select you choice for hot/dry model: ",
                                                c("", "access1_0","bcc_csm1_1","bnu_esm","cccma_canesm2","cesm1_cam5_1_fv2","cnrm_cm5","csiro_mk3_6_0","fio_esm","gfdl_cm3","gfdl_esm2g","gfdl_esm2m","giss_e2_h","giss_e2_r","HadGEM2_AO","inm_cm4","ipsl_cm5a_lr","ipsl_cm5a_mr","lasg_fgoals_g2","miroc_esm","miroc_esm_chem","miroc_miroc5","mohc_hadgem2_cc","mohc_hadgem2_es","mpi_esm_lr","mpi_esm_mr","mri_cgcm3","ncar_ccsm4","ncc_noresm1_m")
                                                , multiple = FALSE, selected = NULL), 
                                      selectInput("cold_wet", "2. Select you choice for cold/wet model: ",
                                                c("", "access1_0","bcc_csm1_1","bnu_esm","cccma_canesm2","cesm1_cam5_1_fv2","cnrm_cm5","csiro_mk3_6_0","fio_esm","gfdl_cm3","gfdl_esm2g","gfdl_esm2m","giss_e2_h","giss_e2_r","HadGEM2_AO","inm_cm4","ipsl_cm5a_lr","ipsl_cm5a_mr","lasg_fgoals_g2","miroc_esm","miroc_esm_chem","miroc_miroc5","mohc_hadgem2_cc","mohc_hadgem2_es","mpi_esm_lr","mpi_esm_mr","mri_cgcm3","ncar_ccsm4","ncc_noresm1_m")
                                                , multiple = FALSE, selected = NULL),
                                      h5("3. Most central model: ", htmlOutput("central", inline=T)),
                                      actionButton("calc.diss", "Calculate dissimilars"),
                                      # h5("4. Most different model to the previous three: ", textOutput("mostdissimilar1", inline=T)),
                                      # h5("5. Most different model to the previous four: ", textOutput("mostdissimilar2", inline=T))#,
                                      h5("4. Most different model to the previous three: ", htmlOutput("mostdissimilar1", inline=T)),
                                      h5("5. Most different model to the previous four: ", htmlOutput("mostdissimilar2", inline=T))#,
                                      ),
                       conditionalPanel(condition = "input.GCMchoice == 'sim'",
                                      h5("These are the models inside a 95% confidence interval"),
                                      dataTableOutput("similars")
                                      )
                     ),
            
            tabPanel("COMPARE - spatial plots",
                     h3("Plots of comparisons with ensembles for each variable"),
                     h5("Violet/positive values indicate a higher prediction than the average
                        among compared GCMs (i.e.: larger increase of temperature/precipitation), yellow/negative
                        values indicate a lower prediction"),
                     h5("(this takes a minute to print, especially if you selected a large area)"),
                     uiOutput("differencePlots"))#,
    ))
  )
))

