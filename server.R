
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
if (!require("pacman")) install.packages("pacman")
pacman::p_load(magrittr, dplyr, raster, ccafs, rasterVis, maptools, shiny, rgeos, sp, viridis, ggplot2, ggforce, rworldxtra, matrixStats)

data(countriesHigh)
# data(wrld_simpl)

# rcp.equiv <- data.frame(name = c("RCP 2.6", "RCP 4.5", "RCP 6.0", "RCP 8.5"), cod = c(26, 45, 60, 85))
# year.equiv <- data.frame(name = c("2050", "2070"), cod = c(50, 70))
my.objects <- reactiveValues(my.extent = extent(-180, 180, -90, 83))
scenarios <- readRDS("scenarios.rds")
guidebioclim <- readRDS("guidebioclim.rds")
my.objects$calculate.dissimilars <- 0

#### Functions
# This function filters only some bioclimatic variables from the set of 19
select.bio <- function(variables, wanted.bio){
  return(lapply(variables, "[[", wanted.bio))
}

# This function crops the variables to the extent of the study area
select.extent <- function(variables, my.extent){
  return(lapply(variables, function(x){crop(x, my.extent)}))
}

# Function that calculates ensembles of chosen models (it will only compare inside scenarios (taken from name))
calculate.ensembleGCM <- function(variables){
  # check which scenarios are being analysed
  my.scenarios <- expand.grid(unique(scenarios$year[scenarios$scenarios %in% names(variables)]), unique(scenarios$rcp[scenarios$scenarios %in% names(variables)]))
  
  my.ensembles <- list()
  for (sc in 1:nrow(my.scenarios)){
    sub.variables <- variables[grepl(as.character(my.scenarios[sc,]$Var1), names(variables))]  # filters year
    sub.variables <- sub.variables[grepl(as.character(my.scenarios[sc,]$Var2), names(sub.variables))] # filters rcp
    if (length(sub.variables) > 3){
      bio.stack <- stack()
      for(b in 1:nlayers(variables[[1]])){
        compare.variables <- stack(lapply(sub.variables, "[[", b))
        ensemble <- mean(compare.variables, na.rm=T) %>% setNames(names(variables[[1]])[b])
        bio.stack <- stack(bio.stack, ensemble)
      }
      my.ensembles[[length(my.ensembles)+1]] <- bio.stack
      names(my.ensembles)[[length(my.ensembles)]] <- paste0("ensemble.", as.character(my.scenarios[sc,]$Var1), ".", as.character(my.scenarios[sc,]$Var2))
    }
  }
  return(my.ensembles)
}

# Compare variables with ensemble
compare.ensembleGCM.scale <- function(variables, ensemble, my.bioclim){#, stat="directional"){
  # Table to store results
  table <- expand.grid(unique(scenarios$year[scenarios$scenarios %in% names(variables)]), unique(scenarios$rcp[scenarios$scenarios %in% names(variables)]), stringsAsFactors = F) %>% arrange(Var1) %>% setNames(c("year", "rcp"))
  index <- data.frame(scenarios[scenarios$scenario %in% names(variables),])
  ensemble.index <- names(ensemble) %>% as.data.frame %>% setNames("ensemble")
  
  clim.names <- my.bioclim
  # if (class(clim.variables) == "character"){
    # if (clim.variables == "all"){clim.names <- names(variables[[1]])}
    # if (clim.variables == "temp"){clim.names <- names(variables[[1]])[names(variables[[1]]) %in% c("bio_1", "bio_2", "bio_3", "bio_4", "bio_5", "bio_6", "bio_7", "bio_8", "bio_9", "bio_10", "bio_11")]}
    # if (clim.variables == "prec"){clim.names <- names(variables[[1]])[names(variables[[1]]) %in% c("bio_12", "bio_13", "bio_14", "bio_15", "bio_16", "bio_17", "bio_18", "bio_19")]}
    # if (clim.variables == c("temp & prec")) {
      # clim.names <- c(names(variables[[1]])[names(variables[[1]]) %in% c("bio_1", "bio_2", "bio_3", "bio_4", "bio_5", "bio_6", "bio_7", "bio_8", "bio_9", "bio_10", "bio_11")], names(variables[[1]])[names(variables[[1]]) %in% c("bio_12", "bio_13", "bio_14", "bio_15", "bio_16", "bio_17", "bio_18", "bio_19")])
    # }
  # }
  # if (class(clim.variables) == "numeric") {clim.names <- paste0("bio_", clim.variables)}
  
  res.list <- list()
  res.list.sd <- list()
  for (s in 1:nrow(table)){
    scen <- table[s,]
    # Ensemble
    scen.ensemble <- ensemble[[paste0("ensemble.", scen$year, ".", scen$rcp)]]
    # Models
    scen.index <- index$scenarios[index$year == scen$year & index$rcp == scen$rcp]
    scen.list <- variables[scen.index]
    
    # Compute result table
    res.table <- index$gcm[index$year == scen$year & index$rcp == scen$rcp] %>% as.data.frame %>% setNames("GCM")
    res.table.sd <- res.table
    # res.table <- as.data.frame(scen.index) %>% setNames("GCM") %>% strsplit(":::")
    for (v in clim.names){
      # Select the bioclim variables for all the models
      v.list <- sapply(scen.list, "[[", v)
      analisis.table <- sapply(v.list, as.vector) %>% t #%>% as.data.frame    # es más rápido como matriz
      # # Mean of absolute values of the GCMs result  -> No se puede hacer con esta aproximación de "scale"
      # if (stat == "absolut"){
      #   analisis.table
      #   table.means <- rowMeans(analisis.table, na.rm=T) %>% as.data.frame %>% setNames(paste0(v)) 
      # }
      # Mean of non-absolute values of the GCMs result, that allows to evaluate directionality
      # if (stat == "directional"){
      table.scale <- scale(analisis.table)
      table.means <- rowMeans(table.scale, na.rm=T) %>% as.data.frame %>% setNames(paste0(v)) 
      res.table <- cbind(res.table, table.means)
      row.names(res.table) <- 1:nrow(res.table)
      
      table.sd <- matrixStats::rowSds(table.scale, na.rm=T) %>% as.data.frame %>% setNames(paste0(v))
      res.table.sd <- cbind(res.table.sd, table.sd)
      row.names(res.table.sd) <- 1:nrow(res.table.sd)
      
      # }
      
    }
    
    res.list[[s]] <- res.table
    names(res.list)[[s]] <- paste0(scen$year, ".", scen$rcp)
    
    res.list.sd[[s]] <- res.table.sd
    names(res.list.sd)[[s]] <- paste0(scen$year, ".", scen$rcp)
  }
  
  # List mean and sd results
  final.res.list <- list()
  final.res.list[[1]] <- res.list
  names(final.res.list)[[1]] <- "Means"
  final.res.list[[2]] <- res.list.sd
  names(final.res.list)[[2]] <- "SD"
  return(final.res.list)
}  

plots.comparisonGCM <- function(variables, ensemble){
  # variables <- clim.vars
  # ensemble <- clim.ensembles
  my.summary <- data.frame(scenario = names(variables), stringsAsFactors = F) %>% 
    cbind(as.data.frame(stringr::str_split(.$scenario, "\\.", simplify=T))) %>% 
    setNames(c("scenarios", "gcm", "year", "rcp")) %>% 
    mutate(year = gsub("year", "", year),
           rcp = gsub("rcp", "", rcp)) 
  
  my.scenarios <- expand.grid(unique(my.summary[,3]), unique(my.summary[,4]))
  
  
  ploted.list <- list()
  
  # Loop to compare each
  for (j in 1:nrow(my.scenarios)){
    sce.gcms <- my.summary[my.summary$year == my.scenarios$Var1[j] & my.summary$rcp == my.scenarios$Var2[j],]
    
    plot.list <- list()
    
    for (k in 1:nrow(sce.gcms)){
      difference <- variables[[sce.gcms$scenarios[k]]] - ensemble[[paste0("ensemble.", sce.gcms$year[k], ".", sce.gcms$rcp[k])]]
      plot.list[[length(plot.list)+1]] <- difference
      names(plot.list)[[length(plot.list)]] <- as.character(sce.gcms$gcm[k])
    }
    # Reorganize layers according to bioclims
    sub.ploted.list <- list()
    for (b in 1:nlayers(difference)){
      b.stack <- stack(lapply(plot.list, function(x){raster::subset(x, b)}))
      # If it is a temperature variable, divide by 10
      if(guidebioclim$transformation[guidebioclim$code == names(plot.list[[1]])[b]] == 1){
        b.stack <- b.stack/10
      }
      
      # Plot
      my.theme <- rasterTheme(pch=19, cex=0.7, region=rev(brewer.pal(9, 'Spectral'))) # azul claro -> rojo
      my.RdBuTheme <- rasterTheme(pch=19, cex=0.7, region= rev(brewer.pal(9, "RdBu")))   # celeste -> rojo
      # res.min2 <- min(values(b.stack), na.rm=T)
      # res.max2 <- max(values(b.stack), na.rm=T)
      res.min2 <- quantile(values(b.stack), .001, na.rm=T)[[1]]
      res.max2 <- quantile(values(b.stack), .999, na.rm=T)[[1]]
      
      highest.res <- which(c(abs(res.min2), abs(res.max2)) == max(abs(res.min2), abs(res.max2)))
      if(highest.res == 1) {res.max2 <- -res.min2}
      if(highest.res == 2) {res.min2 <- -res.max2}
      
      res.at2 <- seq(res.min2, res.max2, length.out=length(my.RdBuTheme$regions$col)-1)
      res.ckey2 <- list(at=res.at2, col=my.RdBuTheme$regions$col)
      bio.comp.plot <- levelplot(b.stack, 
                                 margin=F, 
                                 at=res.at2,
                                 main = paste0(my.scenarios$Var1[j], ", rcp", my.scenarios$Var2[j], ", ", names(plot.list[[1]])[b], " (", guidebioclim$meaning[guidebioclim$code == names(plot.list[[1]])[b]], ")"),
                                 col.regions = c("grey20", rev(plasma(255))),
                                 par.strip.text = list(col = "white"),
                                 par.settings = list(
                                   strip.background = list(col = "grey40")),
                                 maxpixels=1e4)
      #                            colorkey=res.ckey2,
      #                            par.settings=my.RdBuTheme,
      #                            # col.regions = c("grey20", plasma(255)),
      #                            par.strip.text = list(col = "white"),
      #                            # par.settings = list(
      #                            #  strip.background = list(col = "grey40")),
      #                            maxpixels=1e4)
      sub.ploted.list[[length(sub.ploted.list)+1]] <- bio.comp.plot
      names(sub.ploted.list)[[length(sub.ploted.list)]] <- paste0("bio ", names(plot.list[[1]])[b], ".", my.scenarios$Var1[j], ".", my.scenarios$Var2[j])
    }
    ploted.list <- c(ploted.list, sub.ploted.list)
  }  
  return(ploted.list)
}

plots.GCMs <- function(variables){
  # variables <- clim.vars
  
  my.summary <- data.frame(scenario = names(variables), stringsAsFactors = F) %>% 
    cbind(as.data.frame(stringr::str_split(.$scenario, "\\.", simplify=T))) %>% 
    setNames(c("scenarios", "gcm", "year", "rcp")) %>% 
    mutate(year = gsub("year", "", year),
           rcp = gsub("rcp", "", rcp)) 
  
  my.scenarios <- expand.grid(unique(my.summary[,3]), unique(my.summary[,4]))
  
  
  ploted.list <- list()
  
  # Loop to compare each
  for (j in 1:nrow(my.scenarios)){
    sce.gcms <- my.summary[my.summary$year == my.scenarios$Var1[j] & my.summary$rcp == my.scenarios$Var2[j],]
    
    plot.list <- list()
    
    for (k in 1:nrow(sce.gcms)){
      sce.plots <- variables[[sce.gcms$scenarios[k]]]
      plot.list[[length(plot.list)+1]] <- sce.plots
      names(plot.list)[[length(plot.list)]] <- as.character(sce.gcms$gcm[k])
    }
    # Reorganize layers according to bioclims
    sub.ploted.list <- list()
    for (b in 1:nlayers(sce.plots)){
      b.stack <- stack(lapply(plot.list, function(x){raster::subset(x, b)}))
      # If it is a temperature variable, divide by 10
      if(guidebioclim$transformation[guidebioclim$code == names(plot.list[[1]])[b]] == 1){
        b.stack <- b.stack/10
      }
      
      # Plot
      my.theme <- rasterTheme(pch=19, cex=0.7, region=rev(brewer.pal(9, 'Spectral'))) # azul claro -> rojo
      my.RdBuTheme <- rasterTheme(pch=19, cex=0.7, region= rev(brewer.pal(9, "RdBu")))   # celeste -> rojo
      res.min2 <- min(values(b.stack), na.rm=T)
      res.max2 <- max(values(b.stack), na.rm=T)
      # res.min2 <- quantile(values(b.stack), .001, na.rm=T)[[1]]
      # res.max2 <- quantile(values(b.stack), .999, na.rm=T)[[1]]
      
      # highest.res <- which(c(abs(res.min2), abs(res.max2)) == max(abs(res.min2), abs(res.max2)))
      # if(highest.res == 1) {res.max2 <- -res.min2}
      # if(highest.res == 2) {res.min2 <- -res.max2}
      
      res.at2 <- seq(res.min2, res.max2, length.out=length(my.RdBuTheme$regions$col)-1)
      res.ckey2 <- list(at=res.at2, col=my.RdBuTheme$regions$col)
      bio.comp.plot <- levelplot(b.stack, 
                                 margin=F, 
                                 at=res.at2,
                                 main = paste0(my.scenarios$Var1[j], ", rcp", my.scenarios$Var2[j], ", ", names(plot.list[[1]])[b], " (", guidebioclim$meaning[guidebioclim$code == names(plot.list[[1]])[b]], ")"),
                                 col.regions = c("grey20", rev(plasma(255))),
                                 par.strip.text = list(col = "white"),
                                 par.settings = list(
                                   strip.background = list(col = "grey40")),
                                 maxpixels=1e4)
      #                            colorkey=res.ckey2,
      #                            par.settings=my.RdBuTheme,
      #                            # col.regions = c("grey20", plasma(255)),
      #                            par.strip.text = list(col = "white"),
      #                            # par.settings = list(
      #                            #  strip.background = list(col = "grey40")),
      #                            maxpixels=1e4)
      sub.ploted.list[[length(sub.ploted.list)+1]] <- bio.comp.plot
      names(sub.ploted.list)[[length(sub.ploted.list)]] <- paste0("bio ", names(plot.list[[1]])[b], ".", my.scenarios$Var1[j], ".", my.scenarios$Var2[j])
    }
    ploted.list <- c(ploted.list, sub.ploted.list)
  }  
  return(ploted.list)
}


shinyServer(function(input, output) {
  
  # PLOT THE STUDY AREA AND UPDATE IT
  output$distPlot <- renderPlot({
    # generate bins based on input$bins from ui.R
    # plot(crop(wrld_simpl, my.objects$my.extent), col="#2ca25f", border="white")
    plot(crop(countriesHigh, my.objects$my.extent), axes = TRUE, col="#2ca25f", border="white")
    box()
  })
  
  # LOAD AND MAKE CALCULATIONS WITH THE VARIABLES, CREATE ENSEMBLES
  observeEvent(input$updateExtent,{
    cat('\n', 'Updating the extent of the world map')
    ## Get extent
    # Type of extent input
    type <- input$type
    country <- input$country
    
    if (type == "con") { if (country == "World") {my.objects$my.extent <- extent(-180, 180, -90, 83)
    } else {my.objects$my.extent <- countriesHigh[countriesHigh$NAME %in% country,]}
    } else if (type == "num") {my.objects$my.extent <- extent(input$minlon, input$maxlon, input$minlat, input$maxlat)
    } else {e <- input$plot_brush;
    if (is.null(e)) {my.objects$my.extent <- extent(countriesHigh)
    } else {my.objects$my.extent <- extent(e$xmin, e$xmax, e$ymin, e$ymax)}}
    cat(' DONE!')
  })
  
  observeEvent(input$go,{    # When the button "go" is actioned, all comands here inside will be run
    # Load selected variables
    clim.vars.files <- c(scenarios$scenarios[scenarios$year %in% input$year & scenarios$rcp %in% input$rcp & scenarios$gcm %in% input$gcm])
          # year.t <- 2070; rcp.t <- 45; gcm.t <- c("cesm1_cam5_1_fv2","gfdl_cm3","mohc_hadgem2_cc","mohc_hadgem2_es","mpi_esm_mr","bcc_csm1_1","access1_0","ncar_ccsm4","lasg_fgoals_g2","gfdl_esm2g","ncc_noresm1_m","mpi_esm_lr","cnrm_cm5","giss_e2_h")
          # clim.vars.files <- c(scenarios$scenarios[scenarios$year %in% year.t & scenarios$rcp %in% rcp.t & scenarios$gcm %in% gcm.t])
    
    # Reading variables
    cat('\n', 'Loading-in variables')
    withProgress(message = 'Loading GCM variables',
                 detail = 'This should not take long', value = 0,{
                 clim.vars <- list()
                  for (s in clim.vars.files){
                    clim.vars[[length(clim.vars)+1]] <- readRDS(paste0("clim/", scenarios$year[scenarios$scenarios == s], "/rcp", scenarios$rcp[scenarios$scenarios == s], "/", s, ".rds"))
                    names(clim.vars)[[length(clim.vars)]] <- s
                  }
                 })
    cat(' DONE!')
                  
    # Crop climatic variables at defined extent
    cat('\n', 'cropping variables')
    withProgress(message = 'Cropping variables to desired extent',
                 detail = 'If you picked many..., maybe go take a coffee...', value = 0,{
                 cropping.extent <- my.objects$my.extent
                        # extent.t <- c(-110, -34, -57, 29)  
                        # cropping.extent <- extent.t
                  clim.vars <- select.extent(clim.vars, cropping.extent)
                 })
    cat(' DONE!')
    
    # Filter selected bioclimatic variables
    cat('\n', 'selecting bioclimatic variables')
    if (input$analysistype == "tVSp"){
      my.objects$bio.vars <- c(1,5,6,8,9,10,11, 12,13,14,16,17,18,19)
            # bio.vars.t <- c(1,5,6,8,9,10,11, 12,13,14,16,17,18,19)
      clim.vars <- select.bio(clim.vars, my.objects$bio.vars)
            # clim.vars <- select.bio(clim.vars, bio.vars.t)
    }
    if (input$analysistype == "custom"){
      my.objects$bio.vars <- input$selected.bio
      clim.vars <- select.bio(clim.vars, my.objects$bio.vars)
    }
    cat(' DONE!')

    # Calculate ensembles
    cat('\n', 'Calculating ensembles')
    withProgress(message = 'Calculating ensembles',
                 detail = 'This is really not taking so long :)', value = 0,{
                 clim.ensembles <- calculate.ensembleGCM(clim.vars)
                })
    cat(' DONE!')
    
    my.summary <- data.frame(scenario = names(clim.vars), stringsAsFactors = F) %>%
      cbind(as.data.frame(stringr::str_split(.$scenario, "\\.", simplify=T))) %>%
      setNames(c("scenarios", "gcm", "year", "rcp")) %>%
      mutate(year = gsub("year", "", year),
             rcp = gsub("rcp", "", rcp))
    
    
    # Store variables in my.objects box, so they can be used in other observeEvents
    my.objects$clim.vars <- clim.vars
    my.objects$clim.ensembles <- clim.ensembles
    my.objects$my.summary <- my.summary
    my.objects$my.scenarios <- expand.grid(unique(my.summary[,3]), unique(my.summary[,4])) %>% setNames(c("year", "rcp"))
    
    
    ### Plots  -- Creo que mejor en el observeEvent separado... para que no lo haga hasta pasar de pestaña
                  # de todas maneras, aún no lo hace
    # cat('\n', 'preparing differencePlots')
    # differencePlot.list <- plots.ensembleGCM(my.objects$clim.vars, my.objects$clim.ensembles)
        # plot.list <- plots.ensembleGCM(clim.vars, clim.ensembles)
    # cat('DONE')
    # my.objects$differencePlot.list
  })
  ###############################################
  
  output$table <- renderDataTable({
    cat('\n', 'Doing the hard calculations')
    withProgress(message = 'Doing the hard calculations',
                 detail = 'Please, be patient, our best men are on it!', value = 0,{
                  my.bioclim <- paste0("bio_", my.objects$bio.vars)
                        # my.bioclim <- paste0("bio_", bio.vars.t)
                  result.table <- compare.ensembleGCM.scale(my.objects$clim.vars, my.objects$clim.ensembles, my.bioclim)
                        # result.table <- compare.ensembleGCM.scale(clim.vars, clim.ensembles, my.bioclim)
                  result.means <- result.table[["Means"]]
                  result.means2 <- list()
                  result.means2[[1]] <- result.means[[1]]$GCM %>% as.data.frame %>% setNames("GCM") %>% mutate(GCM = as.character(GCM))
                  names(result.means2)[[1]] <- "Temp"
                  for(c in 1:length(result.means)) {result.means2$Temp <- left_join(result.means2$Temp, result.means[[c]][,c(1, which(names(result.means[[c]]) %in% paste0("bio_", 1:11)))], by="GCM")}
                  result.means2[[2]] <- result.means[[1]]$GCM %>% as.data.frame %>% setNames("GCM") %>% mutate(GCM = as.character(GCM))
                  names(result.means2)[[2]] <- "Prec"
                  for(c in 1:length(result.means)) {result.means2$Prec <- left_join(result.means2$Prec, result.means[[c]][,c(1, which(names(result.means[[c]]) %in% paste0("bio_", 12:19)))], by="GCM")}
                  
                  result.means3 <- result.means2[[1]][,1] %>% as.data.frame(stringsAsFactors=F) %>% 
                    cbind(rowMeans(result.means2$Temp[,2:ncol(result.means2$Temp)], na.rm=T)) %>%
                    cbind(rowMeans(result.means2$Prec[,2:ncol(result.means2$Prec)], na.rm=T)) %>% 
                    setNames(c("GCM", "Temp", "Prec")) %>% 
                    mutate(dist = pointDistance(c(0, 0), .[,2:3], lonlat=FALSE)) %>% 
                    arrange(dist)
                  
                  ### SD
                  result.SD <- result.table[["SD"]]
                  result.SD2 <- list()
                  result.SD2[[1]] <- result.SD[[1]]$GCM %>% as.data.frame %>% setNames("GCM") %>% mutate(GCM = as.character(GCM))
                  names(result.SD2)[[1]] <- "Temp"
                  for(c in 1:length(result.SD)) {result.SD2$Temp <- left_join(result.SD2$Temp, result.SD[[c]][,c(1, which(names(result.SD[[c]]) %in% paste0("bio_", 1:11)))], by="GCM")}
                  result.SD2[[2]] <- result.SD[[1]]$GCM %>% as.data.frame %>% setNames("GCM") %>% mutate(GCM = as.character(GCM))
                  names(result.SD2)[[2]] <- "Prec"
                  for(c in 1:length(result.SD)) {result.SD2$Prec <- left_join(result.SD2$Prec, result.SD[[c]][,c(1, which(names(result.SD[[c]]) %in% paste0("bio_", 12:19)))], by="GCM")}
                  
                  result.SD3 <- result.SD2[[1]][,1] %>% as.data.frame(stringsAsFactors=F) %>% 
                    cbind(rowMeans(result.SD2$Temp[,2:ncol(result.SD2$Temp)], na.rm=T)) %>%
                    cbind(rowMeans(result.SD2$Prec[,2:ncol(result.SD2$Prec)], na.rm=T)) %>% 
                    setNames(c("GCM", "Temp", "Prec")) %>% 
                    mutate(dist = pointDistance(c(0, 0), .[,2:3], lonlat=FALSE)) %>% 
                    arrange(dist)
                  })
    cat(' DONE!')
    my.objects$result.means3 <- result.means3
    my.objects$result.SD3 <- result.SD3
    # input$make.differencePlot <- 1 
    
    return(result.means3)
  })
  
  output$pVSt.plot <- renderPlot({
    cat('\n', 'Ploting pVSt')
    scale.max <- max(abs(my.objects$result.means3$Prec), abs(my.objects$result.means3$Temp))+0.2
        # scale.max <- max(abs(result.means3$Prec), abs(result.means3$Temp))+0.2
    
    ## Without confidence interval
    # pVSt.plot <- ggplot(data=my.objects$result.means3, aes(Temp, Prec, label = as.character(GCM))) +
    #     # pVSt.plot <- ggplot(data=result.means3, aes(Temp, Prec, label = as.character(GCM))) +
    #                 # ggtitle("Precipitation Vs Temperature plot") +
    #                 geom_vline(xintercept = 0, color="grey") +
    #                 geom_hline(yintercept = 0, color="grey") +
    #                 geom_point(color="dark blue") +
    #                 coord_fixed() +
    #                 theme_bw() +
    #                 geom_text(size=3, nudge_x = 0.02, nudge_y = 0.02) +
    #                 scale_x_continuous(limits = c(-scale.max,scale.max)) +
    #                 scale_y_continuous(limits = c(-scale.max,scale.max))
    #                 # scale_x_continuous(limits = c((-max(c(max(abs(my.objects$result.means3$Temp)), max(abs(my.objects$result.means3$Prec))))-0.3),(max(c(max(abs(my.objects$result.means3$Temp)), max(abs(my.objects$result.means3$Prec))))+0.3)))
    #                 # scale_y_continuous(limits = c((-max(c(max(abs(my.objects$result.means3$Temp)), max(abs(my.objects$result.means3$Prec))))-0.3),(max(c(max(abs(my.objects$result.means3$Temp)), max(abs(my.objects$result.means3$Prec))))+0.3)))
    
    #### Calculate Confidence Interval
    mean <- mean(my.objects$result.means3$dist, na.rm=T)
    sd <- sd(my.objects$result.means3$dist, na.rm=T)
    serror <- sd/sqrt(nrow(my.objects$result.means3))
    conf.int <- qnorm(.975)*serror
    mean.conf.int <- mean+conf.int
    result.means.sel <- my.objects$result.means3
    result.means.sel$selectionCircle <- ifelse(result.means.sel$dist <= mean.conf.int, 1, 0)
    result.means.sel$selectionCircle <- as.factor(result.means.sel$selectionCircle)
    # plot con confidence
    circle <- data.frame(Temp=0, Prec=0, r=mean.conf.int)

    pVSt.plot <- ggplot() + 
                  ggforce::geom_circle(data=circle, aes(x0=Temp, y0=Prec, r=r, fill=NULL), color="grey") +  # Hay que hacer primero el círculo... si no se lía
                  # ggforce::geom_circle(data=circle2, aes(x0=Temp, y0=Prec, r=r, fill=NULL), color="grey") +  # Hay que hacer primero el círculo... si no se lía
                  # ggtitle(region) +
                  geom_vline(xintercept = 0, color="grey") +
                  geom_hline(yintercept = 0, color="grey") +
                  geom_point(data=result.means.sel, aes(Temp, Prec, label = as.character(GCM), color=selectionCircle)) +#, color="dark blue") +
                  coord_fixed() +
                  # theme_bw() +
                  geom_text(data=result.means.sel, aes(Temp, Prec, label = as.character(GCM)), size=3, nudge_x = 0.1, nudge_y = 0.05) +
                  scale_x_continuous(limits = c(-scale.max,scale.max)) +
                  scale_y_continuous(limits = c(-scale.max,scale.max))
                    
    
    my.objects$start.differencePlot <- "start"
    my.objects$result.means.sel <- result.means.sel
    cat('DONE')
    return(pVSt.plot)
  })
  
  #########################
  ### Selection of models according to the two methods
  output$central <- renderUI({
  # output$central <- renderText({
    # my.objects$result.means3$GCM[my.objects$result.means3$dist == min(my.objects$result.means3$dist, na.rm=T)]
    pplot <- paste("</b>", my.objects$result.means3$GCM[1],"</b>")
    HTML(pplot)
    my.objects$central <- my.objects$result.means3$GCM[1] %>% as.character
          # central <- result.means3$GCM[1] %>% as.character
  })
  
  observeEvent(input$calc.diss, {
  # output$mostdissimilar1 <- renderText({
    hot_dry <- input$hot_dry     #hot_dry <- "cesm1_cam5_1_fv2"
    cold_wet <- input$cold_wet   #cold_wet <- "gfdl_esm2g"
    central <- my.objects$central
    selected123 <- c(hot_dry, cold_wet, central)
    ## Find the forth model (the most different to the other three)
    result.means4 <- my.objects$result.means3
          # result.means4 <- result.means3
    result.means4 <- result.means4 %>% 
        mutate(dist4a = pointDistance(result.means4[result.means4$GCM == selected123[1], 2:3], result.means4[,2:3], lonlat=FALSE)) %>%
        mutate(dist4b = pointDistance(result.means4[result.means4$GCM == selected123[2], 2:3], result.means4[,2:3], lonlat=FALSE)) %>% 
        mutate(dist4c = pointDistance(result.means4[result.means4$GCM == selected123[3], 2:3], result.means4[,2:3], lonlat=FALSE)) %>% 
        mutate(dist4sum = dist4a+dist4b+dist4c) %>% 
        arrange(dist4sum)
    mostdissimilar1 <- result.means4$GCM[nrow(result.means4)]
    if (mostdissimilar1 %in% selected123){mostdissimilar1 <- result.means4$GCM[nrow(result.means4)-1]}
    if (mostdissimilar1 %in% selected123){mostdissimilar1 <- result.means4$GCM[nrow(result.means4)-2]}
    if (mostdissimilar1 %in% selected123){mostdissimilar1 <- result.means4$GCM[nrow(result.means4)-3]}
    # return(mostdissimilar1 %>% as.character)
    selected1234 <- c(selected123, mostdissimilar1)# %>% as.character)
    my.objects$mostdissimilar1 <- mostdissimilar1
  # })#,
  
  # observeEvent(my.objects$selected1234, {
  # output$mostdissimilar2 <- renderText({
    ## Find the fifth model (the most different to the other four)
    # hot_dry <- input$hot_dry
    # cold_wet <- input$cold_wet
    # central <- my.objects$central
    # mostdissimilar1 <- my.objects$mostdissimilar1
    # selected1234 <- c(hot_dry, cold_wet, central, mostdissimilar1 %>% as.character)
    # selected1234 <- my.objects$selected1234

    result.means5 <- my.objects$result.means3
          # result.means5 <- result.means3
    result.means5 <- result.means5 %>%
        mutate(dist5a = pointDistance(result.means5[result.means5$GCM == selected1234[1], 2:3], result.means5[,2:3], lonlat=FALSE)) %>%
        mutate(dist5b = pointDistance(result.means5[result.means5$GCM == selected1234[2], 2:3], result.means5[,2:3], lonlat=FALSE)) %>%
        mutate(dist5c = pointDistance(result.means5[result.means5$GCM == selected1234[3], 2:3], result.means5[,2:3], lonlat=FALSE)) %>%
        mutate(dist5d = pointDistance(result.means5[result.means5$GCM == selected1234[4], 2:3], result.means5[,2:3], lonlat=FALSE)) %>%
        mutate(dist5sum = dist5a+dist5b+dist5c+dist5d) %>%
        arrange(dist5sum)
    mostdissimilar2 <- result.means5$GCM[nrow(result.means5)]
    if (mostdissimilar2 %in% selected1234){mostdissimilar2 <- result.means5$GCM[nrow(result.means5)-1]}
    if (mostdissimilar2 %in% selected1234){mostdissimilar2 <- result.means5$GCM[nrow(result.means5)-2]}
    if (mostdissimilar2 %in% selected1234){mostdissimilar2 <- result.means5$GCM[nrow(result.means5)-3]}
    if (mostdissimilar2 %in% selected1234){mostdissimilar2 <- result.means5$GCM[nrow(result.means5)-4]}
    # return(mostdissimilar2 %>% as.character)
    # my.objects$selected12345 <- c(selected1234, mostdissimilar2)# %>% as.character)
    my.objects$mostdissimilar2 <- mostdissimilar2
  })
  
  # output$mostdissimilar1 <- renderText({
  output$mostdissimilar1 <- renderUI({
    # print(my.objects$selected12345[4])
    that.plot <- paste("<b>",my.objects$mostdissimilar1,"</b>")
    HTML(that.plot)
  })  
  
  # output$mostdissimilar2 <- renderText({
  output$mostdissimilar2 <- renderUI({
    # print(my.objects$selected12345[5])
    this.plot <- paste("<b>",my.objects$mostdissimilar2,"</b>")
    HTML(this.plot)
  })
  
  ## Similar models - 95% confidence interval
  output$similars <- renderDataTable({
    result.means.sel <- my.objects$result.means.sel   
    selected <- result.means.sel$GCM[result.means.sel$selectionCircle == 1]
    selected <- as.data.frame(selected)
    return(selected)
    # this.plot <- paste("<b>",selected,"</b>")
    # HTML(this.plot)
  })
  
  
  
  
  ###############################################
  observeEvent(my.objects$start.differencePlot, {
  # output$differencePlots <- renderPlot({
    cat('\n', 'preparing differencePlots')
    # clim.vars <- my.objects$clim.vars
    # clim.ensembles <- my.objects$clim.ensembles
    withProgress(message = 'Pre-preparing some of the plots',
                 detail = '', value = 0,{
    my.bioclim <- paste0("bio_", my.objects$bio.vars)
    differencePlot.list <- plots.comparisonGCM(my.objects$clim.vars, my.objects$clim.ensembles)
        # differencePlot.list <- plots.comparisonGCM(clim.vars, clim.ensembles)
                 })
       
     ### THIS WILL CREATE AN OUTPUT OBJECT WITH EACH OF THE PLOTS
      for (i in 1:length(differencePlot.list)) {
        # Need local so that each item gets its own number. Without it, the value
        # of i in the renderPlot() will be the same across all instances, because
        # of when the expression is evaluated.
        local({
                my_i <- i
                plotname <- paste0("plot_comp", my_i)
                output[[plotname]] <- renderPlot({
                                          differencePlot.list[[my_i]]
                                          })
              })
      }
    cat('DONE')
    my.objects$differencePlot.list <- differencePlot.list
    # return(differencePlot.list)
  })

  output$differencePlots <- renderUI({
    
    # differencePlot.list <- my.objects$differencePlot.list
    plot_expl_output_list <- lapply(1:length(my.objects$differencePlot.list), function(i) {
    # plot_expl_output_list <- lapply(1:length(differencePlot.list), function(i) {
      plotname <- paste0("plot_comp", i)
      plotOutput(plotname)#, height = 280, width = 1500)
    })
    # Convert the list to a tagList - this is necessary for the list of items
    # to display properly.
    do.call(tagList, plot_expl_output_list)
  })
  
  ###############################################
  observeEvent(input$go, {
  # output$differencePlots <- renderPlot({
    cat('\n', 'preparing GCM Plots')
    # clim.vars <- my.objects$clim.vars
    # clim.ensembles <- my.objects$clim.ensembles
    withProgress(message = 'Pre-preparing some of the plots',
                 detail = '', value = 0,{
    GCMPlot.list <- plots.GCMs(my.objects$clim.vars)
        # GCMPlot.list <- plots.GCMs(clim.vars)
                 })    
     ### THIS WILL CREATE AN OUTPUT OBJECT WITH EACH OF THE PLOTS
      for (i in 1:length(GCMPlot.list)) {
        # Need local so that each item gets its own number. Without it, the value
        # of i in the renderPlot() will be the same across all instances, because
        # of when the expression is evaluated.
        local({
                my_i <- i
                plotname <- paste0("plot_GCM", my_i)
                output[[plotname]] <- renderPlot({
                                          GCMPlot.list[[my_i]]
                                          })
              })
      }
    cat('DONE')
    my.objects$GCMPlot.list <- GCMPlot.list
    # return(GCMPlot.list)
  })

  output$GCMPlots <- renderUI({
    # differencePlot.list <- my.objects$differencePlot.list
    GCMplot_expl_output_list <- lapply(1:length(my.objects$GCMPlot.list), function(i) {
    # GCMplot_expl_output_list <- lapply(1:length(GCMPlot.list), function(i) {
      plotname <- paste0("plot_GCM", i)
      plotOutput(plotname)#, height = 280, width = 1500)
    })
    # Convert the list to a tagList - this is necessary for the list of items
    # to display properly.
    do.call(tagList, GCMplot_expl_output_list)
  })
  ############################################### 
  
  ## Pick models
  
  
  
  
  
  
  
  
  
  # output$visFun <- renderText({c(input$year, input$rcp, input$gcm)})
  # output$info <- renderText({
  #   xy_range_str <- function(e) {
  #     if(is.null(e)) return("NULL\n")
  #     paste0("xmin=", round(e$xmin, 1), " xmax=", round(e$xmax, 1),
  #            " ymin=", round(e$ymin, 1), " ymax=", round(e$ymax, 1))
  #   }
  # 
  #   paste0(
  #     "brush: ", xy_range_str(input$plot_brush)
  #   )
  # 
  # })

  output$Possible <- renderDataTable({
    Possible <- readRDS("Possible.rds")
    as.data.frame(Possible)
  })
  
  # output$table <- renderDataTable({
  # output$vars <- reactive({
  # observeEvent(input$go, {
  #   # vars <- list()
  #   # for (a in input$year){
  #   #   wc.vars.temp3 <- list()
  #   #   for (x in input$rcp){
  #   #     wc.vars.temp2 <- list()
  #   #     withProgress(message = 'Download in progress',
  #   #                  detail = 'This may take a while...', value = 0,
  #   #                  {for (m in input$gcm){
  #   #                    #rm(wc.vars.temp)
  #   #                    # wc.vars.temp <- try(getData("CMIP5", var="bio", res=10, rcp=rcp.equiv$cod[rcp.equiv$name == x], year=year.equiv$cod[year.equiv$name == a], model=m))
  #   #                    wc.vars.temp <- try(getData("CMIP5", var="bio", res=10, rcp=rcp.equiv$cod[rcp.equiv$name == x], year=year.equiv$cod[year.equiv$name == a], model=m, download=T, path="/Volumes/Javier_nuevo/Biologia/Analisis/_Climatic_layers"))
  #   #                    #if (!is.null(wc.vars.temp)){
  #   #                    if (!inherits(wc.vars.temp, "try-error")){     # Some of the combinations are not available, don't stop if you try to download an unexisting combination
  #   #                      names(wc.vars.temp) <- paste0("bio_", 1:19)
  #   #                      wc.vars.temp2[[length(wc.vars.temp2)+1]] <- wc.vars.temp
  #   #                      names(wc.vars.temp2)[[length(wc.vars.temp2)]] <- m
  #   #                    }
  #   #                    #}
  #   #                  }
  #   #                  })
  #   #     wc.vars.temp3[[length(wc.vars.temp3)+1]] <- wc.vars.temp2
  #   #     names(wc.vars.temp3)[[length(wc.vars.temp3)]] <- x
  #   #   }
  #   #   vars[[length(vars)+1]] <- wc.vars.temp3
  #   #   names(vars)[[length(vars)]] <- a
  #   # }
  #   
  #   
  #   ### 2) Filter the variables to a set desired by the user
  #   
  #   ## a) Select only some of the bioclimatic variables
  #   
  #   # Filter selected bios across all the layers
  #   # vars <- vars
  #   # for(y in 1:length(vars)){   # year
  #   #   for (r in 1:length(vars[[y]])){   # rcp
  #   #     for (g in 1:length(vars[[y]][[r]])){   # gcm
  #   #       vars[[y]][[r]][[g]] <- subset(vars[[y]][[r]][[g]], subset=paste0("bio_", input$selected.bio))
  #   #     }
  #   #   }
  #   # }
  #   
  #   
  #   
  #   ## b) Crop variables to a user-defined Extent
  #   
  #   ## cut the variables to the extent if a user-extent is available 
  #   ## (this could be done after bioclim variables have been selected. It would be faster, but not all variables would be ready for an eventual selection later
  #   # vars <- vars
  #   for(y in 1:length(vars)){
  #     for (r in 1:length(vars[[y]])){
  #       for (g in 1:length(vars[[y]][[r]])){
  #         vars[[y]][[r]][[g]] <- if (class(my.objects$my.extent) == "Extent")
  #         {crop(vars[[y]][[r]][[g]], my.objects$my.extent)} else {mask(vars[[y]][[r]][[g]], my.objects$my.extent)}
  #       }
  #     }
  #   }
  #   my.extent$vars <- vars
  # })  
  
    ###############################################
    
    ### 3) Calculate Ensembles and compare each model with them
    ## a) Calculate ensembles
    # output$ensembles <- reactive({  
    # observeEvent(input$go, {
    #   ensembles <- my.objects$clim.vars    # Ensembles will be stored in this object
    #   
    #   for (y in 1:length(my.objects$clim.vars)){        # year
    #     for (r in 1:length(my.objects$clim.vars[[y]])){        # rcp
    #       bio.sub.ensemble <- stack()
    #       for (b in 1:length(input$selected.bio)){         # for every bio-variable
    #         gcm.sub.ensemble <- stack()
    #         for (g in 1:length(my.objects$clim.vars[[y]][[r]])){      # across gcms
    #           gcm.sub.ensemble <- stack(gcm.sub.ensemble, my.objects$clim.vars[[y]][[r]][[g]][[b]])   
    #         }
    #         # Calculate the mean across gcms
    #         bio.sub.ensemble <- stack(bio.sub.ensemble, mean(gcm.sub.ensemble, na.rm = TRUE) %>% setNames(paste0("bio_", input$selected.bio[b])))
    #       }
    #       ensembles[[y]][[r]] <- bio.sub.ensemble    # and store it in "ensembles"
    #     }
    #   }
    #   my.extent$ensembles <- ensembles
    # })
    
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
  
    # output$table <- renderDataTable({
    #   comp.table <- data.frame(year = character(), rcp = character(), gcm = character())
    #   comp.table$year <- as.character(comp.table$year); comp.table$rcp <- as.character(comp.table$rcp); comp.table$gcm <- as.character(comp.table$gcm)
    #   for (b in input$selected.bio){
    #     comp.table$newcol <- numeric(nrow(comp.table))
    #     names(comp.table)[ncol(comp.table)] <- paste0("bio_", b)
    #   }
    #   comp.table.template <- comp.table.net <- comp.table.abs <- comp.table
    #   
    #   # Loop throw all bioclims of each GCM in each scenario and compare it to the correspondent ensemble
    #   
    #   plot.list <- my.objects$clim.vars
    #   
    #   for(y in 1:length(my.extent$ensembles)){
    #     for (r in 1:length(my.extent$ensembles[[y]])){
    #       for (g in 1:length(my.extent$vars[[y]][[r]])){
    #         comp.table.temp <- comp.table.template
    #         comp.table.temp[nrow(comp.table.temp)+1,1] <- names(my.extent$ensembles)[[y]]    # these lines prepare the data for the gcm info to bind to the table
    #         comp.table.temp[nrow(comp.table.temp),2] <- names(my.extent$ensembles[[y]])[[r]]
    #         comp.table.temp[nrow(comp.table.temp),3] <- names(my.extent$vars[[y]][[r]])[[g]]
    #         
    #         # Raster stack with difference rasters between this GCM layers and ensembles'
    #         difference <- my.extent$vars[[y]][[r]][[g]] - my.extent$ensembles[[y]][[r]]
    #         
    #         # Sum values (absolute and net)
    #         res.abs <- cellStats(abs(difference), stat="sum", na.rm=T)
    #         res.net <- cellStats(difference, stat="sum", na.rm=T)
    #         
    #         # Save to tables
    #         comp.table.temp.abs <- comp.table.temp.net <- comp.table.temp
    #         
    #         comp.table.temp.abs[nrow(comp.table.temp.abs), 4:(3+length(res.abs))] <- res.abs
    #         comp.table.abs <- rbind(comp.table.abs, comp.table.temp.abs)
    #         
    #         comp.table.temp.net[nrow(comp.table.temp.net), 4:(3+length(res.net))] <- res.net
    #         comp.table.net <- rbind(comp.table.net, comp.table.temp.net)
    #         
    #         ### Save difference plot in the list to provide the final plot
    #         plot.list[[y]][[r]][[g]] <- difference
    #       }
    #     }
    #   }  
    #   
    #   ### Normalize the results, so different variables can be compared among them
    #   # a) Tabla de valores absolutos
    #   normalizeMinMax <- function(x, newMin, newMax){ (x - min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T)) * (newMax - newMin) + newMin }
    #   comp.table.abs.norm <- comp.table.abs
    #   
    #   for (y in unique(comp.table.abs.norm$year)){
    #     for (r in unique(comp.table.abs.norm$rcp)){
    #       for (b in grep("bio", names(comp.table.abs.norm))){
    #         sc.values <- comp.table.abs.norm[(comp.table.abs.norm$year == y & comp.table.abs.norm$rcp == r),b]
    #         norm.values <- normalizeMinMax(sc.values, 0, 1)
    #         comp.table.abs.norm[(comp.table.abs.norm$year == y & comp.table.abs.norm$rcp == r),b] <- norm.values
    #       }
    #     }
    #   }
    #   # Summarize result in a TOTAL final column
    #   for (row in 1:nrow(comp.table.abs.norm)){
    #     comp.table.abs.norm$total[row] <- mean(as.numeric(comp.table.abs.norm[row, grep("bio", names(comp.table.abs.norm))]), na.rm=T)
    #   }
    #   # comp.table.abs.norm
    #   
    #   
    #   # b) Tabla de valores netos
    #   normalizeMinMax_Net <- function(x, newMin, newMax){
    #     (x - dataMin)/(dataMax-dataMin) * (newMax - newMin) + newMin }
    #   
    #   comp.table.net.norm <- comp.table.net
    #   
    #   for (y in unique(comp.table.net.norm$year)){
    #     for (r in unique(comp.table.net.norm$rcp)){
    #       for (b in grep("bio", names(comp.table.net.norm))){
    #         sc.values <- comp.table.net.norm[(comp.table.net.norm$year == y & comp.table.net.norm$rcp == r),b]
    #         
    #         highest <- which(c(abs(max(sc.values)), abs(min(sc.values))) == max(c(abs(max(sc.values)), abs(min(sc.values))))) # = 1 if positive values are greater, =2 if negative are greater  
    #         if(highest == 1) {dataMax <- max(sc.values); dataMin <- -max(sc.values)}
    #         if(highest == 2) {dataMin <- min(sc.values); dataMax <- -min(sc.values)}
    #         
    #         norm.values <- normalizeMinMax_Net(sc.values, -1, 1)
    #         comp.table.net.norm[(comp.table.net.norm$year == y & comp.table.net.norm$rcp == r),b] <- norm.values
    #       }
    #     }
    #   }
    #   # Summarize result in a TOTAL final column
    #   for (row in 1:nrow(comp.table.net.norm)){
    #     comp.table.net.norm$total[row] <- mean(as.numeric(comp.table.net.norm[row, grep("bio", names(comp.table.net.norm))]), na.rm=T)
    #   }
    #   # print(comp.table.net.norm)
    #   comp.table.norm <- comp.table.net.norm
    #   return(comp.table.norm)
    # })  
      

    # ### PLOTS!
    # 
    # # Plot chosen variables
    # observeEvent(input$go, {
    # # output$explore.plot <- renderPlot({
    #   explore.plot <- my.extent$vars
    #   list.expl.plots <- list()
    #   
    #   for (y in names(explore.plot)){
    #     for (r in names(explore.plot[[y]])){
    #       names.list <- names(explore.plot[[y]][[r]])
    #       # bio.plot.stack <- vector("list", length(names.list))
    #       for (b in 1:length(input$selected.bio)){
    #         b.stack <- stack()
    #         for (g in names(explore.plot[[y]][[r]])){
    #           g.layer <- explore.plot[[y]][[r]][[g]][[b]] %>% setNames(g)
    #           b.stack <- stack(b.stack, g.layer)
    #         }  
    #         ## Plot it!
    #         
    #         my.theme <- rasterTheme(pch=19, cex=0.7, region=rev(brewer.pal(9, 'Spectral'))) # azul claro -> rojo
    #         my.RdBuTheme <- rasterTheme(pch=19, cex=0.7, region= rev(brewer.pal(9, "RdBu")))   # celeste -> rojo
    #         res.min <- min(values(b.stack), na.rm=T)
    #         res.max <- max(values(b.stack), na.rm=T)
    #         # Si hay outliers, puede verse mal... por eso antes usaba esto del cuartil
    #         # res.min <- quantile(values(b.stack), .001, na.rm=T)[[1]]
    #         # res.max <- quantile(values(b.stack), .999, na.rm=T)[[1]]
    #         
    #         res.at <- seq(res.min, res.max, length.out=length(my.RdBuTheme$regions$col)-1)
    #         # res.ckey <- list(at=res.at, col=my.RdBuTheme$regions$col)
    #         
    #         b.plot <- levelplot(b.stack, 
    #                             margin=F, 
    #                             at=res.at, 
    #                             main = paste0("GCM versions of bio ", input$selected.bio[b]),
    #                             col.regions = c("grey20", viridis(255)),
    #                             par.strip.text = list(col = "white"),
    #                             par.settings = list(
    #                               strip.background = list(col = "grey40")),
    #                             maxpixels=1e4)
    #         
    #         # bio.plot.stack[[b]] <- b.plot
    #         list.expl.plots[[length(list.expl.plots)+1]] <- b.plot
    #       }
    #       # names(bio.plot.stack) <- paste0("bio_", input$selected.bio)
    #       # explore.plot[[y]][[r]] <- bio.plot.stack
    #     }
    #   }
    #   # explore.plot
    #   my.extent$list.expl.plots <- list.expl.plots
    #   
    #   ### THIS WILL CREATE AN OUTPUT OBJECT WITH EACH OF THE PLOTS
    #   for (i in 1:length(list.expl.plots)) {
    #     # Need local so that each item gets its own number. Without it, the value
    #     # of i in the renderPlot() will be the same across all instances, because
    #     # of when the expression is evaluated.
    #     local({
    #       my_i <- i
    #       plotname <- paste0("plot_expl", my_i)
    #       output[[plotname]] <- renderPlot({
    #         list.expl.plots[[my_i]]
    #       })
    #     })
    #   }
    # })
    # 
    # 
    # 
    # observeEvent(input$go, {
    # # output$compare.plot <- renderPlot({
    #   # plot.list.2 <- plot.list <- my.extent$vars
    #   compare.plot <- my.extent$vars
    #   list.of.plots <- list()
    #   for (y in names(compare.plot)){
    #     for (r in names(compare.plot[[y]])){
    #       names.list <- names(compare.plot[[y]][[r]])
    #       # bio.compare.plot.stack <- vector("list", length(names.list))
    #       for (b in 1:length(input$selected.bio)){
    #         bio.comp.stack <- stack()
    #         for (g in names(compare.plot[[y]][[r]])){
    #           difference <- compare.plot[[y]][[r]][[g]][[b]] - my.extent$ensembles[[y]][[r]][[b]]
    #           bio.comp.stack <- stack(bio.comp.stack, difference %>% setNames(g))
    #         }  
    #         ## Plot it!
    #         
    #         my.theme <- rasterTheme(pch=19, cex=0.7, region=rev(brewer.pal(9, 'Spectral'))) # azul claro -> rojo
    #         my.RdBuTheme <- rasterTheme(pch=19, cex=0.7, region= rev(brewer.pal(9, "RdBu")))   # celeste -> rojo
    #         res.min2 <- min(values(bio.comp.stack), na.rm=T)
    #         res.max2 <- max(values(bio.comp.stack), na.rm=T)
    #         # res.min2 <- quantile(values(bio.comp.stack), .001, na.rm=T)[[1]]
    #         # res.max2 <- quantile(values(bio.comp.stack), .999, na.rm=T)[[1]]
    #         
    #         highest.res <- which(c(abs(res.min2), abs(res.max2)) == max(abs(res.min2), abs(res.max2)))
    #         if(highest.res == 1) {res.max2 <- -res.min2}
    #         if(highest.res == 2) {res.min2 <- -res.max2}
    #         
    #         res.at2 <- seq(res.min2, res.max2, length.out=length(my.RdBuTheme$regions$col)-1)
    #         res.ckey2 <- list(at=res.at2, col=my.RdBuTheme$regions$col)
    #         
    #         bio.comp.plot <- levelplot(bio.comp.stack, 
    #                                    margin=F, 
    #                                    at=res.at2,
    #                                    colorkey=res.ckey2,
    #                                    main = paste0(y, ", ", r, ", bio ", input$selected.bio[b]),
    #                                    par.settings=my.RdBuTheme,
    #                                    # col.regions = c("grey20", plasma(255)),
    #                                    par.strip.text = list(col = "white"),
    #                                    # par.settings = list(
    #                                    #  strip.background = list(col = "grey40")),
    #                                    maxpixels=1e4)
    #         
    #         list.of.plots[[length(list.of.plots)+1]] <- bio.comp.plot
    #       }
    #       # names(bio.compare.plot.stack) <- paste0("bio_", input$selected.bio)
    #       # compare.plot[[y]][[r]] <- bio.compare.plot.stack
    #     }
    #   }
    #   my.extent$list.of.plots <- list.of.plots
    #   
    #   ### THIS WILL CREATE AN OUTPUT OBJECT WITH EACH OF THE PLOTS
    #   for (i in 1:length(list.of.plots)) {
    #     # Need local so that each item gets its own number. Without it, the value
    #     # of i in the renderPlot() will be the same across all instances, because
    #     # of when the expression is evaluated.
    #     local({
    #       my_i <- i
    #       plotname <- paste0("plot", my_i)
    #       output[[plotname]] <- renderPlot({
    #         list.of.plots[[my_i]]
    #       })
    #     })
    #   }
    # })
    # 
    # 
    # ### THIS USES renderUI TO PRINT AS MANY PLOTS AS THERE ARE (***explore.plot)
    # output$explore.plot <- renderUI({
    #   plot_expl_output_list <- lapply(1:length(my.extent$list.expl.plots), function(i) {
    #     plotname <- paste0("plot_expl", i)
    #     plotOutput(plotname)#, height = 280, width = 1500)
    #   })
    #   # Convert the list to a tagList - this is necessary for the list of items
    #   # to display properly.
    #   do.call(tagList, plot_expl_output_list)
    # })
    # 
    # 
    # ### THIS USES renderUI TO PRINT AS MANY PLOTS AS THERE ARE (***compare.plot)
    # output$compare.plot <- renderUI({
    #   plot_output_list <- lapply(1:length(my.extent$list.of.plots), function(i) {
    #     plotname <- paste0("plot", i)
    #     plotOutput(plotname)#, height = 280, width = 1500)
    #   })
    #   # Convert the list to a tagList - this is necessary for the list of items
    #   # to display properly.
    #   do.call(tagList, plot_output_list)
    # })
    
    
  
  
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



