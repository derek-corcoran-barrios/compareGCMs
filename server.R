
# This is the server logic for a Shiny web application.
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
data(wrld_simpl)
rcp.equiv <- data.frame(name = c("RCP 2.6", "RCP 4.5", "RCP 6.0", "RCP 8.5"), cod = c(26, 45, 60, 85))
year.equiv <- data.frame(name = c("2050", "2070"), cod = c(50, 70))
my.extent <- reactiveValues(EXT = extent(-180, 180, -90, 83))

shinyServer(function(input, output) {

   observeEvent(input$go,{
    type <- input$type
    country <- input$country
    
    if (type == "con") { if (country == "World") {my.extent$EXT <- extent(-180, 180, -90, 83)
    } else {my.extent$EXT <- wrld_simpl[wrld_simpl$NAME %in% country,]}
    } else if (type == "num") {my.extent$EXT <- extent(input$minlon, input$maxlon, input$minlat, input$maxlat)
    } else {e <- input$plot_brush;
    if (is.null(e)) {my.extent$EXT <- extent(wrld_simpl)
    } else {my.extent$EXT <- extent(e$xmin, e$xmax, e$ymin, e$ymax)}}
  })
  
  # PLOT THE STUDY AREA AND UPDATE IT
  output$distPlot <- renderPlot({
    # generate bins based on input$bins from ui.R
    plot(crop(wrld_simpl, my.extent$EXT))
    
  })
  
  output$visFun <- renderText({c(input$year, input$rcp, input$all.models)})
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
  
  output$Possible <- renderDataTable({
    Possible <- readRDS("Possible.rds")
    as.data.frame(Possible)
  })
  
  output$table <- renderDataTable({
    vars <- list()
    for (a in input$year){
      wc.vars.temp3 <- list()
      for (x in input$rcp){
        wc.vars.temp2 <- list()
        withProgress(message = 'Download in progress',
                     detail = 'This may take a while...', value = 0,
        {for (m in input$all.models){
          #rm(wc.vars.temp)
          wc.vars.temp <- try(getData("CMIP5", var="bio", res=10, rcp=rcp.equiv$cod[rcp.equiv$name == x], year=year.equiv$cod[year.equiv$name == a], model=m))
          #if (!is.null(wc.vars.temp)){
            if (!inherits(wc.vars.temp, "try-error")){     # Some of the combinations are not available, don't stop if you try to download an unexisting combination
              names(wc.vars.temp) <- paste0("bio_", 1:19)
              wc.vars.temp2[[length(wc.vars.temp2)+1]] <- wc.vars.temp
              names(wc.vars.temp2)[[length(wc.vars.temp2)]] <- m
            }
          #}
        }
        })
        wc.vars.temp3[[length(wc.vars.temp3)+1]] <- wc.vars.temp2
        names(wc.vars.temp3)[[length(wc.vars.temp3)]] <- x
      }
      vars[[length(vars)+1]] <- wc.vars.temp3
      names(vars)[[length(vars)]] <- a
    }
    
    ### 2) Filter the variables to a set desired by the user
    
    ## a) Select only some of the bioclimatic variables
    
    # Filter selected bios across all the layers
    vars <- vars
    for(y in 1:length(vars)){   # year
      for (r in 1:length(vars[[y]])){   # rcp
        for (g in 1:length(vars[[y]][[r]])){   # gcm
          vars[[y]][[r]][[g]] <- subset(vars[[y]][[r]][[g]], subset=paste0("bio_", input$selected.bio))
        }
      }
    }
    
    
    
    ## b) Crop variables to a user-defined Extent
    
    ## cut the variables to the extent if a user-extent is available 
    ## (this could be done after bioclim variables have been selected. It would be faster, but not all variables would be ready for an eventual selection later
      vars <- vars
      for(y in 1:length(vars)){
        for (r in 1:length(vars[[y]])){
          for (g in 1:length(vars[[y]][[r]])){
            vars[[y]][[r]][[g]] <- if (class(my.extent$EXT) == "Extent")
              {crop(vars[[y]][[r]][[g]], my.extent$EXT)} else {trim(mask(vars[[y]][[r]][[g]], my.extent$EXT))}
          }
        }
      }
    
    
    ###############################################
    
    ### 3) Calculate Ensembles and compare each model with them
    ## a) Calculate ensembles
    
    ensembles <- vars    # Ensembles will be stored in this object
    
    for (y in 1:length(vars)){        # year
      for (r in 1:length(vars[[y]])){        # rcp
        bio.sub.ensemble <- stack()
        for (b in 1:length(input$selected.bio)){         # for every bio-variable
          gcm.sub.ensemble <- stack()
          for (g in 1:length(vars[[y]][[r]])){      # across gcms
            gcm.sub.ensemble <- stack(gcm.sub.ensemble, vars[[y]][[r]][[g]][[b]])   
          }
          # Calculate the mean across gcms
          bio.sub.ensemble <- stack(bio.sub.ensemble, mean(gcm.sub.ensemble, na.rm = TRUE) %>% setNames(paste0("bio_", input$selected.bio[b])))
        }
        ensembles[[y]][[r]] <- bio.sub.ensemble    # and store it in "ensembles"
      }
    }
    
    ## b) Compare each GCM variable with the ensemble
    # Create a table to store all the comparisons
    
    comp.table <- data.frame(year = character(), rcp = character(), gcm = character())
    comp.table$year <- as.character(comp.table$year); comp.table$rcp <- as.character(comp.table$rcp); comp.table$gcm <- as.character(comp.table$gcm)
    for (b in input$selected.bio){
      comp.table$newcol <- numeric(nrow(comp.table))
      names(comp.table)[ncol(comp.table)] <- paste0("bio_", b)
    }
    comp.table.template <- comp.table
    
    # Loop throw all bioclims of each GCM in each scenario and compare it to the correspondent ensemble
    for(y in 1:length(ensembles)){
      for (r in 1:length(ensembles[[y]])){
        for (g in 1:length(vars[[y]][[r]])){
          comp.table.temp <- comp.table.template
          comp.table.temp[nrow(comp.table.temp)+1,1] <- names(ensembles)[[y]]    # these lines prepare the data for the gcm info to bind to the table
          comp.table.temp[nrow(comp.table.temp),2] <- names(ensembles[[y]])[[r]]
          comp.table.temp[nrow(comp.table.temp),3] <- names(vars[[y]][[r]])[[g]]
          res <- cellStats(abs(vars[[y]][[r]][[g]] - ensembles[[y]][[r]]), stat="sum", na.rm=TRUE)  # Calculate the sum of the differences (in absolute value)
          comp.table.temp[nrow(comp.table.temp), 4:(3+length(res))] <- res
          comp.table <- rbind(comp.table, comp.table.temp)
        }
      }
    }  
    
    ### Normalize the results, so different variables can be compared among them
    # Function
    normalizeMinMax <- function(x, newMin, newMax){ (x - min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T)) * (newMax - newMin) + newMin }
    
    comp.table.norm <- comp.table
    
    for (y in unique(comp.table.norm$year)){
      for (r in unique(comp.table.norm$rcp)){
        for (b in grep("bio", names(comp.table.norm))){
          sc.values <- comp.table.norm[(comp.table.norm$year == y & comp.table.norm$rcp == r),b]
          norm.values <- normalizeMinMax(sc.values, 0, 1)
          comp.table.norm[(comp.table.norm$year == y & comp.table.norm$rcp == r),b] <- norm.values
        }
      }
    }
    
    
    ### Summarize the result: The final TOTAL column summarizes how similar is, in average, the layer to the ensemble
    # The closest to 0 it is, the more similar
    for (row in 1:nrow(comp.table.norm)){
      comp.table.norm$total[row] <- mean(as.numeric(comp.table.norm[row, grep("bio", names(comp.table.norm))]), na.rm=TRUE)
    }
    as.data.frame(comp.table.norm)
  })
  
  output$myYeartabs = renderUI({
    myTabs = lapply(paste('Year', input$year) , tabPanel)
    do.call(tabsetPanel, myTabs)
  })
  
  output$myBiotabs = renderUI({
    myTabs = lapply(paste('Bio', input$selected.bio), tabPanel)
    do.call(tabsetPanel, myTabs)
  })
  
})
