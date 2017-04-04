# Libraries, wd and data

if (!require("pacman")) install.packages("pacman")

pacman::p_load(magrittr, dplyr, dtplyr, raster, openxlsx, ccafs, rasterVis, maptools, shiny)



####### The analysis is composed by 4 parts:
# 1) Read the variables that will be analyzed on a list of raster stacks
# 2) Filter and crop the variables to a set desired by the user (only a few bioclims, extent of your study - area)
# 3) Calculate Ensembles and compare each model with them


####################
# First, set some PARAMETERS of the models you want to compare and the study area
year <- c("2050", "2070")    # from a list with 2050, 2070
rcp <- c("RCP 4.5", "RCP 8.5")   # from a list with "RCP 2.6", "RCP 4.5", "RCP 6.0", "RCP 8.5"

# USER EXTENT
# a. Define the extent directly
my.extent <- extent(-120, -31, -57, 33)

# b. Choose it with drawExtent
plot(vars[[1]][[1]][[1]][[1]]) # Make a plot of one of the variables to help the drawing
my.extent <- drawExtent()

# c. Take it from a country name
country <- "Argentina"
if (class(my.extent) == "character"){
  data(wrld_simpl)
  my.extent <- wrld_simpl[wrld_simpl$NAME %in% country,] %>%
    extent
}
####################


### 1) Read the variables that will be analyzed on a list of raster stacks

# There are 3 options: 
# a) download variables form worldclim
# b) download variables form CCAFS (more complite, sometimes the server is down)
# c) Make you own stacks of variables (see at the end of this section how to organize them on a list)

### a) Download from WorldClim
rcp.equiv <- data.frame(name = c("RCP 2.6", "RCP 4.5", "RCP 6.0", "RCP 8.5"), cod = c(26, 45, 60, 85))
year.equiv <- data.frame(name = c("2050", "2070"), cod = c(50, 70))

#all.models <- c("AC", "BC", "CC", "CE", "CN", "GF", "GD", "GS", "HD", "HG", "HE", "IN", "IP", "MI", "MR", "MC", "MP", "MG", "NO")

all.models <- c("AC", "BC")

year <- c("2050")    # from a list with 2050, 2070
rcp <- c("RCP 4.5", "RCP 8.5")


vars <- list()
for (a in year){
  wc.vars.temp3 <- list()
  for (x in rcp){
    wc.vars.temp2 <- list()
    for (m in all.models){
      rm(wc.vars.temp)
      wc.vars.temp <- try(getData("CMIP5", var="bio", res=10, rcp=rcp.equiv$cod[rcp.equiv$name == x], year=year.equiv$cod[year.equiv$name == a], model=m))
      if (!is.null(wc.vars.temp)){
        if (!inherits(wc.vars.temp, "try-error")){     # Some of the combinations are not available, don't stop if you try to download an unexisting combination
          names(wc.vars.temp) <- paste0("bio_", 1:19)
          wc.vars.temp2[[length(wc.vars.temp2)+1]] <- wc.vars.temp
          names(wc.vars.temp2)[[length(wc.vars.temp2)]] <- m
        }
      }
    }
    wc.vars.temp3[[length(wc.vars.temp3)+1]] <- wc.vars.temp2
    names(wc.vars.temp3)[[length(wc.vars.temp3)]] <- x
  }
  vars[[length(vars)+1]] <- wc.vars.temp3
  names(vars)[[length(wc.vars)]] <- a
}



### b) Download variable rasters from CCAFS
# Table of correspondence between parameters and CCAFS ids
rcp.equiv <- data.frame(name = c("RCP 2.6", "RCP 4.5", "RCP 6.0", "RCP 8.5"), cod = 7:10)
year.equiv <- data.frame(name = c("2050", "2070"), cod = c(6, 8))

# Find URLs with download links (Parameter options with ?'ccafs-search')
cc.query <- list()
for (a in 1:length(year)){
  for (x in 1:length(rcp)){
    cc.query.temp <- cc_search(file_set = 12,
                               scenario = rcp.equiv$cod[rcp.equiv$name == rcp[x]],    # 7:rcp26 ; 8:rcp45 ; 9:rcp60 ; 10:rcp85
                               # model = 2,     # Here you can choose along GCMs, blank will return all available
                               extent = "global", 
                               format = "ascii", 
                               period = year.equiv$cod[year.equiv$name == rcp[a]],       #5:2040 ; 6:2050 ; 7:2060 ; 8:2070
                               variable = 1,     # 1:bioclimatic
                               resolution = 4    # 4:10min
    )# %>% print
    cc.query[length(cc.query)+1] <- list(cc.query.temp)
    names(cc.query)[[x]] <- rcp[x]
  }
}

# Download the variables in those URL and read them as raster stacks
# Store all them in the list "vars"
vars <- list()
for (x in 1:length(cc.query)){
  vars.temp <- lapply(cc.query[[x]][1:length(cc.query[[x]])], cc_data_fetch, progress = TRUE) %>%
    lapply(cc_data_read)
  
  # Change the index of each element in the list by the GCM name
  splitted.url <- strsplit(cc.query[[x]], "/")
  for (u in 1:length(splitted.url)){
    long <- splitted.url[u][[1]] %>%
      tail(n=1)
    short <- strsplit(long, "_rcp")[[1]] %>%
      head(n=1)
    names(vars.temp)[[u]] <- short
  }
  
  vars[[x]] <- vars.temp
  names(vars)[[x]] <- rcp[x]
}   


### c) Alternativelly, the variables could be read from raster files locally.
## For the rest of the script to work straightforward, keep in mind that they should be organized as follows:
## A list, with an element for each year
## Each year in the list, should contain an element for each RCP
## Each RCP should contain stacks of each GCM bioclimatic variables (e.g: a stack for bnu_esm containing bio_1:bio_19), year and rcp given by its position on the list

###############################################

### 2) Filter the variables to a set desired by the user

## a) Select only some of the bioclimatic variables
selected.bio <- c(1:5, 12:13)

# Filter selected bios across all the layers
vars.2 <- vars
for(y in 1:length(vars.2)){   # year
  for (r in 1:length(vars.2[[y]])){   # rcp
    for (g in 1:length(vars.2[[y]][[r]])){   # gcm
      vars.2[[y]][[r]][[g]] <- subset(vars.2[[y]][[r]][[g]], subset=paste0("bio_", selected.bio))
    }
  }
}



## b) Crop variables to a user-defined Extent

## cut the variables to the extent if a user-extent is available 
## (this could be done after bioclim variables have been selected. It would be faster, but not all variables would be ready for an eventual selection later
if (exists("my.extent")){
  vars.3 <- vars.2
  for(y in 1:length(vars.3)){
    for (r in 1:length(vars.3[[y]])){
      for (g in 1:length(vars.3[[y]][[r]])){
        vars.3[[y]][[r]][[g]] <- crop(vars.3[[y]][[r]][[g]], my.extent)
      }
    }
  }
}


###############################################

### 3) Calculate Ensembles and compare each model with them
## a) Calculate ensembles

ensembles <- vars.3    # Ensembles will be stored in this object

for (y in 1:length(vars.3)){        # year
  for (r in 1:length(vars.3[[y]])){        # rcp
    bio.sub.ensemble <- stack()
    for (b in 1:length(selected.bio)){         # for every bio-variable
      gcm.sub.ensemble <- stack()
      for (g in 1:length(vars.3[[y]][[r]])){      # across gcms
        gcm.sub.ensemble <- stack(gcm.sub.ensemble, vars.3[[y]][[r]][[g]][[b]])   
      }
      # Calculate the mean across gcms
      bio.sub.ensemble <- stack(bio.sub.ensemble, mean(gcm.sub.ensemble) %>% setNames(paste0("bio_", selected.bio[b])))
    }
    ensembles[[y]][[r]] <- bio.sub.ensemble    # and store it in "ensembles"
  }
}

## b) Compare each GCM variable with the ensemble
# Create a table to store all the comparisons

comp.table <- data.frame(year = character(), rcp = character(), gcm = character())
comp.table$year <- as.character(comp.table$year); comp.table$rcp <- as.character(comp.table$rcp); comp.table$gcm <- as.character(comp.table$gcm)
for (b in selected.bio){
  comp.table$newcol <- numeric(nrow(comp.table))
  names(comp.table)[ncol(comp.table)] <- paste0("bio_", b)
}
comp.table.template <- comp.table

# Loop throw all bioclims of each GCM in each scenario and compare it to the correspondent ensemble
for(y in 1:length(ensembles)){
  for (r in 1:length(ensembles[[y]])){
    for (g in 1:length(vars.3[[y]][[r]])){
      comp.table.temp <- comp.table.template
      comp.table.temp[nrow(comp.table.temp)+1,1] <- names(ensembles)[[y]]    # these lines prepare the data for the gcm info to bind to the table
      comp.table.temp[nrow(comp.table.temp),2] <- names(ensembles[[y]])[[r]]
      comp.table.temp[nrow(comp.table.temp),3] <- names(vars.3[[y]][[r]])[[g]]
      res <- cellStats(abs(vars.3[[y]][[r]][[g]] - ensembles[[y]][[r]]), stat="sum", na.rm=TRUE)  # Calculate the sum of the differences (in absolute value)
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
  comp.table.norm$total[row] <- mean(as.numeric(comp.table.norm[row, grep("bio", names(comp.table.norm))]), na.rm=T)
}
print(comp.table.norm)


(0.5*5)+(0.3*5)+(1.2*0.2)

