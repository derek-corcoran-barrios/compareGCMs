
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
library(viridis)
library(rworldxtra)
data(countriesHigh)
rcp.equiv <- data.frame(name = c("RCP 2.6", "RCP 4.5", "RCP 6.0", "RCP 8.5"), cod = c(26, 45, 60, 85))
year.equiv <- data.frame(name = c("2050", "2070"), cod = c(50, 70))
my.extent <- reactiveValues(EXT = extent(-180, 180, -90, 83))
my.vars <- reactiveValues(vars = NULL, ensembles = NULL, list.of.plots = NULL, list.expl.plots = NULL)

shinyServer(function(input, output) {
  
  observeEvent(input$go,{
    type <- input$type
    country <- input$country
    
    if (type == "con") { if (country == "World") {my.extent$EXT <- extent(-180, 180, -90, 83)
    } else {my.extent$EXT <- countriesHigh[countriesHigh$ADMIN %in% country,]}
    } else if (type == "num") {my.extent$EXT <- extent(input$minlon, input$maxlon, input$minlat, input$maxlat)
    } else {e <- input$plot_brush;
    if (is.null(e)) {my.extent$EXT <- extent(countriesHigh)
    } else {my.extent$EXT <- extent(e$xmin, e$xmax, e$ymin, e$ymax)}}
  })
  
  # PLOT THE STUDY AREA AND UPDATE IT
  output$distPlot <- renderPlot({
    # generate bins based on input$bins from ui.R
    plot(crop(countriesHigh, my.extent$EXT), axes = TRUE)
    box()
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
  
  # output$table <- renderDataTable({
  # output$vars <- reactive({
  observeEvent(input$go, {
    vars <- list()
    for (a in input$year){
      wc.vars.temp3 <- list()
      for (x in input$rcp){
        wc.vars.temp2 <- list()
        withProgress(message = 'Download in progress',
                     detail = 'This may take a while, 1/12', value = 0,
                     {for (m in input$all.models){
                       #rm(wc.vars.temp)
                       wc.vars.temp <- try(getData("CMIP5", var="bio", res=10, rcp=rcp.equiv$cod[rcp.equiv$name == x], year=year.equiv$cod[year.equiv$name == a], model=m))
                       # wc.vars.temp <- try(getData("CMIP5", var="bio", res=10, rcp=rcp.equiv$cod[rcp.equiv$name == x], year=year.equiv$cod[year.equiv$name == a], model=m, download=T, path="/Volumes/Javier_nuevo/Biologia/Analisis/_Climatic_layers"))
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
    # vars <- vars
    withProgress(message = 'Processing',
                 detail = 'This may take a while, 2/12', value = 1,
                 {
    for(y in 1:length(vars)){   # year
      for (r in 1:length(vars[[y]])){   # rcp
        for (g in 1:length(vars[[y]][[r]])){   # gcm
          vars[[y]][[r]][[g]] <- subset(vars[[y]][[r]][[g]], subset=paste0("bio_", input$selected.bio))
        }
      }
    }})
    
    
    
    ## b) Crop variables to a user-defined Extent
    
    ## cut the variables to the extent if a user-extent is available 
    ## (this could be done after bioclim variables have been selected. It would be faster, but not all variables would be ready for an eventual selection later
    # vars <- vars
    withProgress(message = 'Processing',
                 detail = 'This may take a while, 3/12', value = 1,
                 {
    
    for(y in 1:length(vars)){
      for (r in 1:length(vars[[y]])){
        for (g in 1:length(vars[[y]][[r]])){
          vars[[y]][[r]][[g]] <- if (class(my.extent$EXT) == "Extent")
          {crop(vars[[y]][[r]][[g]], my.extent$EXT)} else {crop(mask(vars[[y]][[r]][[g]], my.extent$EXT),extent(my.extent$EXT))}
        }
      }
    }})
    my.extent$vars <- vars
  })  
  
    ###############################################
    
    ### 3) Calculate Ensembles and compare each model with them
    ## a) Calculate ensembles
    # output$ensembles <- reactive({  
    observeEvent(input$go, {
      ensembles <- my.extent$vars    # Ensembles will be stored in this object
      withProgress(message = 'Processing',
                   detail = 'This may take a while, 4/12', value = 1,
                   {
      for (y in 1:length(my.extent$vars)){        # year
        for (r in 1:length(my.extent$vars[[y]])){        # rcp
          bio.sub.ensemble <- stack()
          for (b in 1:length(input$selected.bio)){         # for every bio-variable
            gcm.sub.ensemble <- stack()
            for (g in 1:length(my.extent$vars[[y]][[r]])){      # across gcms
              gcm.sub.ensemble <- stack(gcm.sub.ensemble, my.extent$vars[[y]][[r]][[g]][[b]])   
            }
            # Calculate the mean across gcms
            bio.sub.ensemble <- stack(bio.sub.ensemble, mean(gcm.sub.ensemble, na.rm = TRUE) %>% setNames(paste0("bio_", input$selected.bio[b])))
          }
          ensembles[[y]][[r]] <- bio.sub.ensemble    # and store it in "ensembles"
        }
      }})
      my.extent$ensembles <- ensembles
    })
    
    # ## b) Compare each GCM variable with the ensemble
    # # Create a table to store all the comparisons
    # 
    # comp.table <- data.frame(year = character(), rcp = character(), gcm = character())
    # comp.table$year <- as.character(comp.table$year); comp.table$rcp <- as.character(comp.table$rcp); comp.table$gcm <- as.character(comp.table$gcm)
    # for (b in input$selected.bio){
    #   comp.table$newcol <- numeric(nrow(comp.table))
    #   names(comp.table)[ncol(comp.table)] <- paste0("bio_", b)
    # }
    # comp.table.template <- comp.table
    # 
    # # Loop throw all bioclims of each GCM in each scenario and compare it to the correspondent ensemble
    # for(y in 1:length(ensembles)){
    #   for (r in 1:length(ensembles[[y]])){
    #     for (g in 1:length(vars[[y]][[r]])){
    #       comp.table.temp <- comp.table.template
    #       comp.table.temp[nrow(comp.table.temp)+1,1] <- names(ensembles)[[y]]    # these lines prepare the data for the gcm info to bind to the table
    #       comp.table.temp[nrow(comp.table.temp),2] <- names(ensembles[[y]])[[r]]
    #       comp.table.temp[nrow(comp.table.temp),3] <- names(vars[[y]][[r]])[[g]]
    #       res <- cellStats(abs(vars[[y]][[r]][[g]] - ensembles[[y]][[r]]), stat="sum", na.rm=TRUE)  # Calculate the sum of the differences (in absolute value)
    #       comp.table.temp[nrow(comp.table.temp), 4:(3+length(res))] <- res
    #       comp.table <- rbind(comp.table, comp.table.temp)
    #     }
    #   }
    # }  
    # 
    # ### Normalize the results, so different variables can be compared among them
    # # Function
    # normalizeMinMax <- function(x, newMin, newMax){ (x - min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T)) * (newMax - newMin) + newMin }
    # 
    # comp.table.norm <- comp.table
    # 
    # for (y in unique(comp.table.norm$year)){
    #   for (r in unique(comp.table.norm$rcp)){
    #     for (b in grep("bio", names(comp.table.norm))){
    #       sc.values <- comp.table.norm[(comp.table.norm$year == y & comp.table.norm$rcp == r),b]
    #       norm.values <- normalizeMinMax(sc.values, 0, 1)
    #       comp.table.norm[(comp.table.norm$year == y & comp.table.norm$rcp == r),b] <- norm.values
    #     }
    #   }
    # }
    # 
    # 
    # ### Summarize the result: The final TOTAL column summarizes how similar is, in average, the layer to the ensemble
    # # The closest to 0 it is, the more similar
    # for (row in 1:nrow(comp.table.norm)){
    #   comp.table.norm$total[row] <- mean(as.numeric(comp.table.norm[row, grep("bio", names(comp.table.norm))]), na.rm=TRUE)
    # }
    # as.data.frame(comp.table.norm)
    
  
    ####################################################################################################
    ## b) Compare each GCM variable with the ensemble
    # Create a table to store all the comparisons
    
  
  
    output$table <- renderDataTable({
      #validate(need(my.vars$vars,"Please, press Update variables"))
      comp.table <- data.frame(year = character(), rcp = character(), gcm = character())
      comp.table$year <- as.character(comp.table$year); comp.table$rcp <- as.character(comp.table$rcp); comp.table$gcm <- as.character(comp.table$gcm)
      for (b in input$selected.bio){
        comp.table$newcol <- numeric(nrow(comp.table))
        names(comp.table)[ncol(comp.table)] <- paste0("bio_", b)
      }
      comp.table.template <- comp.table.net <- comp.table.abs <- comp.table
      
      # Loop throw all bioclims of each GCM in each scenario and compare it to the correspondent ensemble
      
      plot.list <- my.extent$vars
      withProgress(message = 'Processing',
                   detail = 'This may take a while, 5/12', value = 1,
                   {
      for(y in 1:length(my.extent$ensembles)){
        for (r in 1:length(my.extent$ensembles[[y]])){
          for (g in 1:length(my.extent$vars[[y]][[r]])){
            comp.table.temp <- comp.table.template
            comp.table.temp[nrow(comp.table.temp)+1,1] <- names(my.extent$ensembles)[[y]]    # these lines prepare the data for the gcm info to bind to the table
            comp.table.temp[nrow(comp.table.temp),2] <- names(my.extent$ensembles[[y]])[[r]]
            comp.table.temp[nrow(comp.table.temp),3] <- names(my.extent$vars[[y]][[r]])[[g]]
            
            # Raster stack with difference rasters between this GCM layers and ensembles'
            difference <- my.extent$vars[[y]][[r]][[g]] - my.extent$ensembles[[y]][[r]]
            
            # Sum values (absolute and net)
            res.abs <- cellStats(abs(difference), stat="sum", na.rm=T)
            res.net <- cellStats(difference, stat="sum", na.rm=T)
            
            # Save to tables
            comp.table.temp.abs <- comp.table.temp.net <- comp.table.temp
            
            comp.table.temp.abs[nrow(comp.table.temp.abs), 4:(3+length(res.abs))] <- res.abs
            comp.table.abs <- rbind(comp.table.abs, comp.table.temp.abs)
            
            comp.table.temp.net[nrow(comp.table.temp.net), 4:(3+length(res.net))] <- res.net
            comp.table.net <- rbind(comp.table.net, comp.table.temp.net)
            
            ### Save difference plot in the list to provide the final plot
            plot.list[[y]][[r]][[g]] <- difference
          }
        }
      } })
      
      ### Normalize the results, so different variables can be compared among them
      # a) Tabla de valores absolutos
      normalizeMinMax <- function(x, newMin, newMax){ (x - min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T)) * (newMax - newMin) + newMin }
      comp.table.abs.norm <- comp.table.abs
      withProgress(message = 'Processing',
                   detail = 'This may take a while, 6/12', value = 1,
                   {
      for (y in unique(comp.table.abs.norm$year)){
        for (r in unique(comp.table.abs.norm$rcp)){
          for (b in grep("bio", names(comp.table.abs.norm))){
            sc.values <- comp.table.abs.norm[(comp.table.abs.norm$year == y & comp.table.abs.norm$rcp == r),b]
            norm.values <- normalizeMinMax(sc.values, 0, 1)
            comp.table.abs.norm[(comp.table.abs.norm$year == y & comp.table.abs.norm$rcp == r),b] <- norm.values
          }
        }
      }})
      # Summarize result in a TOTAL final column
      withProgress(message = 'Processing',
                   detail = 'This may take a while, 7/12', value = 1,
                   {
      for (row in 1:nrow(comp.table.abs.norm)){
        comp.table.abs.norm$total[row] <- mean(as.numeric(comp.table.abs.norm[row, grep("bio", names(comp.table.abs.norm))]), na.rm=T)
      }})
      # comp.table.abs.norm
      
      
      # b) Tabla de valores netos
      normalizeMinMax_Net <- function(x, newMin, newMax){
        (x - dataMin)/(dataMax-dataMin) * (newMax - newMin) + newMin }
      
      comp.table.net.norm <- comp.table.net
      withProgress(message = 'Processing',
                   detail = 'This may take a while, 8/12', value = 1,
                   {
      for (y in unique(comp.table.net.norm$year)){
        for (r in unique(comp.table.net.norm$rcp)){
          for (b in grep("bio", names(comp.table.net.norm))){
            sc.values <- comp.table.net.norm[(comp.table.net.norm$year == y & comp.table.net.norm$rcp == r),b]
            
            highest <- which(c(abs(max(sc.values)), abs(min(sc.values))) == max(c(abs(max(sc.values)), abs(min(sc.values))))) # = 1 if positive values are greater, =2 if negative are greater  
            if(highest == 1) {dataMax <- max(sc.values); dataMin <- -max(sc.values)}
            if(highest == 2) {dataMin <- min(sc.values); dataMax <- -min(sc.values)}
            
            norm.values <- normalizeMinMax_Net(sc.values, -1, 1)
            comp.table.net.norm[(comp.table.net.norm$year == y & comp.table.net.norm$rcp == r),b] <- norm.values
          }
        }
      }})
      # Summarize result in a TOTAL final column
      for (row in 1:nrow(comp.table.net.norm)){
        comp.table.net.norm$total[row] <- mean(as.numeric(comp.table.net.norm[row, grep("bio", names(comp.table.net.norm))]), na.rm=T)
      }
      # print(comp.table.net.norm)
      comp.table.norm <- comp.table.net.norm
      comp.table.norm$AbsValue <- comp.table.abs.norm$total 
      return(comp.table.norm)
    })  
      
  
    ### PLOTS!
    
    # Plot chosen variables
    observeEvent(input$go, {
    # output$explore.plot <- renderPlot({
      explore.plot <- my.extent$vars
      list.expl.plots <- list()
      withProgress(message = 'Processing',
                   detail = 'This may take a while, 9/12', value = 1,
                   {
      for (y in names(explore.plot)){
        for (r in names(explore.plot[[y]])){
          names.list <- names(explore.plot[[y]][[r]])
          # bio.plot.stack <- vector("list", length(names.list))
          for (b in 1:length(input$selected.bio)){
            b.stack <- stack()
            for (g in names(explore.plot[[y]][[r]])){
              g.layer <- explore.plot[[y]][[r]][[g]][[b]] %>% setNames(g)
              b.stack <- stack(b.stack, g.layer)
            }  
            ## Plot it!
            
            my.theme <- rasterTheme(pch=19, cex=0.7, region=rev(brewer.pal(9, 'Spectral'))) # azul claro -> rojo
            my.RdBuTheme <- rasterTheme(pch=19, cex=0.7, region= rev(brewer.pal(9, "RdBu")))   # celeste -> rojo
            res.min <- min(values(b.stack), na.rm=T)
            res.max <- max(values(b.stack), na.rm=T)
            # Si hay outliers, puede verse mal... por eso antes usaba esto del cuartil
            # res.min <- quantile(values(b.stack), .001, na.rm=T)[[1]]
            # res.max <- quantile(values(b.stack), .999, na.rm=T)[[1]]
            
            res.at <- seq(res.min, res.max, length.out=length(my.RdBuTheme$regions$col)-1)
            # res.ckey <- list(at=res.at, col=my.RdBuTheme$regions$col)
            
            b.plot <- levelplot(b.stack, 
                                margin=F, 
                                at=res.at, 
                                # main = paste0("GCM versions of bio ", input$selected.bio[b]),
                                main = paste0(y, ", ", r, ", bio ", input$selected.bio[b]),
                                col.regions = c("grey20", viridis(255)),
                                par.strip.text = list(col = "white"),
                                par.settings = list(
                                  strip.background = list(col = "grey40")),
                                maxpixels=1e4)
            
            # bio.plot.stack[[b]] <- b.plot
            list.expl.plots[[length(list.expl.plots)+1]] <- b.plot
          }
          # names(bio.plot.stack) <- paste0("bio_", input$selected.bio)
          # explore.plot[[y]][[r]] <- bio.plot.stack
        }
      }})
      # explore.plot
      my.extent$list.expl.plots <- list.expl.plots
      
      ### THIS WILL CREATE AN OUTPUT OBJECT WITH EACH OF THE PLOTS
      withProgress(message = 'Ploting...',
                   detail = 'This may take a while, 10/12', value = 1,
                   {
      for (i in 1:length(list.expl.plots)) {
        # Need local so that each item gets its own number. Without it, the value
        # of i in the renderPlot() will be the same across all instances, because
        # of when the expression is evaluated.
        local({
          my_i <- i
          plotname <- paste0("plot_expl", my_i)
          output[[plotname]] <- renderPlot({
            list.expl.plots[[my_i]]
          })
        })
      }})
    })
    
    
    
    observeEvent(input$go, {
    # output$compare.plot <- renderPlot({
      # plot.list.2 <- plot.list <- my.extent$vars
      compare.plot <- my.extent$vars
      list.of.plots <- list()
      withProgress(message = 'Ploting',
                   detail = 'This may take a while, 11/12', value = 1,
                   {
      for (y in names(compare.plot)){
        for (r in names(compare.plot[[y]])){
          names.list <- names(compare.plot[[y]][[r]])
          # bio.compare.plot.stack <- vector("list", length(names.list))
          for (b in 1:length(input$selected.bio)){
            bio.comp.stack <- stack()
            for (g in names(compare.plot[[y]][[r]])){
              difference <- compare.plot[[y]][[r]][[g]][[b]] - my.extent$ensembles[[y]][[r]][[b]]
              bio.comp.stack <- stack(bio.comp.stack, difference %>% setNames(g))
            }  
            ## Plot it!
            
            my.theme <- rasterTheme(pch=19, cex=0.7, region=rev(brewer.pal(9, 'Spectral'))) # azul claro -> rojo
            my.RdBuTheme <- rasterTheme(pch=19, cex=0.7, region= rev(brewer.pal(9, "RdBu")))   # celeste -> rojo
            res.min2 <- min(values(bio.comp.stack), na.rm=T)
            res.max2 <- max(values(bio.comp.stack), na.rm=T)
            # res.min2 <- quantile(values(bio.comp.stack), .001, na.rm=T)[[1]]
            # res.max2 <- quantile(values(bio.comp.stack), .999, na.rm=T)[[1]]
            
            highest.res <- which(c(abs(res.min2), abs(res.max2)) == max(abs(res.min2), abs(res.max2)))
            if(highest.res == 1) {res.max2 <- -res.min2}
            if(highest.res == 2) {res.min2 <- -res.max2}
            
            res.at2 <- seq(res.min2, res.max2, length.out=length(my.RdBuTheme$regions$col)-1)
            res.ckey2 <- list(at=res.at2, col=my.RdBuTheme$regions$col)
            
            bio.comp.plot <- levelplot(bio.comp.stack, 
                                       margin=F, 
                                       at=res.at2,
                                       colorkey=res.ckey2,
                                       main = paste0(y, ", ", r, ", bio ", input$selected.bio[b]),
                                       par.settings=my.RdBuTheme,
                                       # col.regions = c("grey20", rev(plasma(255))),
                                       par.strip.text = list(col = "white"),
                                       # par.settings = list(
                                       #  strip.background = list(col = "grey40")),
                                       maxpixels=1e4)
            
            list.of.plots[[length(list.of.plots)+1]] <- bio.comp.plot
          }
          # names(bio.compare.plot.stack) <- paste0("bio_", input$selected.bio)
          # compare.plot[[y]][[r]] <- bio.compare.plot.stack
        }
      }})
      my.extent$list.of.plots <- list.of.plots
      
      ### THIS WILL CREATE AN OUTPUT OBJECT WITH EACH OF THE PLOTS
      for (i in 1:length(list.of.plots)) {
        # Need local so that each item gets its own number. Without it, the value
        # of i in the renderPlot() will be the same across all instances, because
        # of when the expression is evaluated.
        local({
          my_i <- i
          plotname <- paste0("plot", my_i)
          output[[plotname]] <- renderPlot({
            list.of.plots[[my_i]]
          })
        })
      }
    })
    
    
    ### THIS USES renderUI TO PRINT AS MANY PLOTS AS THERE ARE (***explore.plot)
    output$explore.plot <- renderUI({
      withProgress(message = 'Ploting',
                   detail = 'This may take a while, 12/12', value = 1,
                   {
      plot_expl_output_list <- lapply(1:length(my.extent$list.expl.plots), function(i) {
        plotname <- paste0("plot_expl", i)
        plotOutput(plotname)#, height = 280, width = 1500)
      })
      # Convert the list to a tagList - this is necessary for the list of items
      # to display properly.
      do.call(tagList, plot_expl_output_list)})
    })
    

    ### THIS USES renderUI TO PRINT AS MANY PLOTS AS THERE ARE (***compare.plot)
    output$compare.plot <- renderUI({
      plot_output_list <- lapply(1:length(my.extent$list.of.plots), function(i) {
        plotname <- paste0("plot", i)
        plotOutput(plotname)#, height = 280, width = 1500)
      })
      # Convert the list to a tagList - this is necessary for the list of items
      # to display properly.
      do.call(tagList, plot_output_list)
    })
    
    
  
  
  #################################
  
  
  # output$myYeartabs = renderUI({
  #   myTabs = lapply(paste('Year', input$year) , tabPanel)
  #   do.call(tabsetPanel, myTabs)
  # })
  # 
  # output$myBiotabs = renderUI({
  #   myTabs = lapply(paste('Bio', input$selected.bio), tabPanel)
  #   do.call(tabsetPanel, myTabs)
  # })
  
})



