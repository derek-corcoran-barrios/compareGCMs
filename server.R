
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyServer(function(input, output) {

  output$distPlot <- renderPlot({
    # generate bins based on input$bins from ui.R
    data(wrld_simpl)
    
    rcp.equiv <- data.frame(name = c("RCP 2.6", "RCP 4.5", "RCP 6.0", "RCP 8.5"), cod = c(26, 45, 60, 85))
    year.equiv <- data.frame(name = c("2050", "2070"), cod = c(50, 70))
    type <- input$type
    country <- input$country
    plot(wrld_simpl)
    e <- input$plot_brush
    
    
    if(type == "click"){
      plot(wrld_simpl)
      if(is.null(e)){plot(wrld_simpl)}
         else{
         my.extent <- extent(e$xmin, e$xmax,e$ymin,e$ymax)
         plot(crop(wrld_simpl, my.extent))}
    } else if(type == "con"){
    
        if(country == "World"){
          plot(wrld_simpl)
        } else {
        my.extent <- wrld_simpl[wrld_simpl$NAME %in% country,] %>% extent
        plot(crop(wrld_simpl, my.extent))
    }} else {
      my.extent <- extent(input$minlon, input$maxlon, input$minlat, input$maxlat)
      plot(crop(wrld_simpl, my.extent))
    }
  })
  
  output$visFun <- renderPrint({c(input$year, input$rcp, input$all.models)})
  output$info <- renderText({
    xy_range_str <- function(e) {
      if(is.null(e)) return("NULL\n")
      paste0("xmin=", round(e$xmin, 1), " xmax=", round(e$xmax, 1), 
             " ymin=", round(e$ymin, 1), " ymax=", round(e$ymax, 1))
    }
    
    paste0(
      "brush: ", xy_range_str(input$plot_brush)
    )
  })
  
})
