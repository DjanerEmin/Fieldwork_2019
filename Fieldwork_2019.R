# Fieldwork 2019
###########################################################################################################################
# Intro
# The HyPlant fieldowrk took place on 26 and 27 June 2019
# Methods
# Experimental setup:
#       One transect with 10 locations
#       Two days of measurements
#       SIF was measured 6 times per location throughout the course of a single day. (the transect was repeated 6 times)
#                   Each measurement lasted 5 min (at location 1 - 20 minutes)
#                   => one transect sampled in 9 x 5 + 20 = 65 minutes
#       LAI was measured 2 times per location - in the morning and in the afternoon.
#       
# Results
#       SIF - in total 10 locations x 6 reps x 2 days = 120 measurements 
#       LAI - in total 10 locations x 2 reps x 2 days = 40 measurements
#       
#       figure 1 LAI estimate changes during the day
#       figure 2 LAI - SIF
#       figure 3 LAI - delta SIF
#       
# Discussion
# 

########################################################################################################################
# Functions

# Extract GPS from FLOX
# this is a function to extract the gps x and y coordinates from the raw FLOX data
# Inputs: dir - the path of raw flox data. It must contain the "date" folders (e.g. 190626, 190627, ...)
# Returns a df containing 
#                       Local and UTC time 
#                       x and y coordinates in WGS 84 Geographical coordinate system (decimal degrees)
extract_gps <- function(dir){ 
    dirs <- list.dirs(path = dir) # specify path to raw or processed flox data
    datalist <- list()
    
    for (i in 2:length(dirs)) {
        # defining pattern is tricky, first take the capital CSV files
        files <- list.files(dirs[i],  pattern = ".CSV")         # ...then use the second half (QE Pro data) 
        files <- files[seq(length(files)/2 + 1, length(files))] # it doesn't really matter 
        #read all those CSV into a list and clean up
        rawlist <- lapply(paste(dirs[i], "/", files, sep = ""), 
                          read.csv, 
                          sep = ";", 
                          header = FALSE)                       # and clean up
        n.list <- lapply(rawlist, function(x) {x<- x[seq(1, nrow(x), 6),c(1,3,17,21,23)]})
        #append to a bigger list for all the folders
        datalist[[i-1]] <- plyr::rbind.fill(n.list)
        
    }
    
    # more cleaning up, remove the N and E from the gps strings
    datalist <- lapply(datalist, function(x) {cbind(x, a <-as.numeric(substr(x$V21, 1, 8)))})
    datalist <- lapply(datalist, function(x) {cbind(x, a <-as.numeric(substr(x$V23, 1, 8)))})
    
    # put some names and substract what is needed
    datalist <- lapply(datalist, function(x) {
        names(x)<- c("V1","localT", "utc", "V2", "V3", "y", "x") 
        x[,c(1,2,3,6,7)]})
    
    # unlist
    nn.df <- plyr::rbind.fill(datalist)
    nn.df$V1 <- as.numeric(nn.df$V1)
    
    return(nn.df)
    #rm(datalist, list, n.list, nn.df, rawlist, dirs, files, i)
    
}


# Center of Gravity
# This is a function to calculate the mean coordinates of a point pattern,
# also known as the center of gravity, such as the one from diurnal FLOX data
# It excludes outliers based on user defined interval "probs"
# Inputs: x1, y1 - Long and Lat respectively 
#         by - grouping (in a newer version this will be optional)
#         probs - the quantile range to consider,e.g. c(.5,.95) means that 5 % of data points will be 
#                   excluded from each end of the distribution
# Returns a data frame

center.of.gravity <- function(x1, y1, by, probs = c(.05,.95)){
    lv <- levels(by)
    a <- list()
    for (i in 1:length(lv)){
        y <- y1[by == lv[i]]
        x <- x1[by == lv[i]]
        qy <- quantile(y, probs = probs)
        qx <- quantile(x, probs = probs)  
        
        a[[i]] <- data.frame(by = lv[i],
                             y = mean(y[y > qy[1] & y < qy[2]]),
                             x= mean(x[x > qx[1] & x < qx[2]]))
        
    }
    
    return(do.call(rbind.data.frame, a))
    
}

#Function to convert Lat/Long (coming from FLOX) to UTM32N used in HyPlant

library(sp)
library(rgdal)

latlong_to_UTM<-function(x,y,zone){
    xy <- data.frame(ID = 1:length(x), X = x, Y = y)
    coordinates(xy) <- c("X", "Y")
    proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
    res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
    return(as.data.frame(res))
}


# RGB image plotter for ggplot

rggbplot <- function(inRGBRst,npix=NA,scale = 'lin'){
    
    rgblinstretch <- function(rgbDf){
        maxList <- apply(rgbDf,2,max)
        minList <- apply(rgbDf,2,min)
        temp<-rgbDf
        for(i in c(1:3)){
            temp[,i] <- (temp[,i]-minList[i])/(maxList[i]-minList[i])
        }
        return(temp)
    }
    
    rgbeqstretch<-function(rgbDf){
        
        temp<-rgbDf
        for(i in c(1:3)){
            unique <- na.omit(temp[,i])
            if (length(unique>0)){
                ecdf<-ecdf(unique)
                temp[,i] <- apply(temp[,i,drop=FALSE],2,FUN=function(x) ecdf(x))
            }
        }
        return(temp)
    }
    
    if(is.na(npix)){
        if(raster::ncell(inRGBRst)>5000){
            npix <- 5000
        }
        else{
            npix <- raster::ncell(inRGBRst)
        }
    }
    x <- raster::sampleRegular(inRGBRst, size=npix, asRaster = TRUE)
    dat <- as.data.frame(as(x, "SpatialPixelsDataFrame"), xy=TRUE)
    colnames(dat)[1:3]<-c('r','g','b')
    
    if(scale=='lin'){
        dat[,1:3]<- rgblinstretch(dat[,1:3])
    } else if(scale=='stretch'){
        dat[,1:3]<- rgbeqstretch(dat[,1:3])
    }
    
    p <- ggplot()+ geom_tile(data=dat, aes(x=x, y=y, fill=rgb(r,g,b))) + scale_fill_identity()
    
}
##########################################################################################################################
# Intro
# The HyPlant fieldowrk took place on 26 and 27 Jue 2019
# Methods
# Experimental setup:
#       One transect with 10 locations
#       Two days of measurements
#       SIF was measured 6 times per location throughout the course of a single day. (the transect was repeated 6 times)
#                   Each measurement lasted 5 min (at location 1 - 20 minutes)
#                   => one transect sampled in 9 x 5 + 20 = 65 minutes
#       LAI was measured 2 times per location - in the morning and in the afternoon.
#       
# Results
#       SIF - in total 10 locations x 6 reps x 2 days = 120 measurements 
#       LAI - in total 10 locations x 2 reps x 2 days = 40 measurements
#       
#       figure 1 LAI estimate changes during the day
#       figure 2 LAI - SIF
#       figure 3 LAI - delta SIF
#       
# Discussion
# 

##########################################################################################################################
### Plotting constants   
pos.nudge <- rnorm(10)/70
size = 2 
stroke = 1.5

##########################################################################################################################
### Input
require(tidyverse)
## load SunScan data and preprocess
suns <- read_delim('sunscan_lai.csv', 
           delim = ',',
           skip = 11) %>% 
    mutate(Day = if_else(row_number() <239,
                         as.Date('2019-06-26'),
                         as.Date('2019-06-27')),
           Rep = if_else(SZA < 45,2,1)) %>%
    group_by(Day, Rep, Location) %>%
    summarise_all(funs(mean, sd)) 


# select LAI only
lai <- suns %>% select(Day, Rep, Location, LAI_mean, LAI_sd)
lai

## test wether the change in LAI from Morning to Afternoon is significant 
# using Paired t-test
lai.tt <- lai %>% select(Location, LAI_mean, Rep) %>% spread(Rep, value = LAI_mean)
t.test(lai.tt$`1`, lai.tt$`2`, alternative = "greater", paired = TRUE)
rm(lai.tt)

## calcualte delta LAI (change between morning and afternoon measurement of LAI)
# for
days <- c('2019-06-26', '2019-06-27')

# initiate a list
deltalai <- list()

for(day in days){
    deltalai[[day]] <- 
        lai %>% 
            select(-LAI_sd) %>%  #for LAI_mean only
            filter(Day == day) %>% 
            group_by(Rep) %>% 
            spread(key = 'Location',value = LAI_mean)
        
}
# calculate delta
deltalai <- lapply(deltalai, function(x) x[2,-c(1,2)] - x[1, -c(1,2)] )
deltalai

deltalai <- as_tibble(t(do.call(rbind, deltalai))) %>%
    gather('date', 'delta') %>%
    bind_cols(lai %>% filter(Rep == 1))
# A tibble: 20 x 2

## plot the LAI change between morning and afternoon Sunscan measurements 

lai2 <- 
    lai %>% 
    mutate(Rep2 = if_else(Rep == 1, '1_Morning', '2_Afternoon'),
               Date = as.factor(Day),
               Location2 = as.character(Location)) 


plotlist <- list()
# plot 26-06
for(day in days){
    plotlist[[paste(day, "LAI")]] <- 
    lai2 %>%
        filter(Day == day) %>%
        ggplot(aes(Rep2, LAI_mean, color = Location2)) +
        geom_point(position =  position_nudge(x = pos.nudge , y = 0)) + 
        geom_segment(aes(x = rep('1_Morning', 10), xend = rep('2_Afternoon', 10), y= `1`, yend = `2`),
                     position =  position_nudge(x = pos.nudge , y = 0), 
                     data = lai %>% 
                         filter( Day == day) %>% 
                         mutate(Location2 = as.character(Location)) %>%
                         select(Location2,LAI_mean, Rep) %>% 
                         spread(Rep, LAI_mean)) + 
        geom_errorbar(aes(ymin = LAI_mean - LAI_sd, ymax = LAI_mean + LAI_sd),
                      width=.005,
                      position =  position_nudge(x = pos.nudge , y = 0))
    
}

rm(lai2, day)


####################################################################################################################
#### FLOX
# mapping the field

selras <- list()
selras$trcol <- raster::stack("icos_color_composites")[[3:1]]
selras$facol <- raster::stack("icos_color_composites")[[4:2]]
# if needed rename

new_names <- substr(names(selras$trcol), 13, 24)
selras <- lapply(selras, setNames, nm = new_names)

# cut to sugar beet extend
win <-
    raster::extent(320600.6, 320800.0, 5638100, 5639550) # the sugar beet field
selras <- lapply(selras, function(x) raster::crop(x, win))



# Map
plotlist[["True Color Map"]] <- rggbplot(selras$trcol, npix = 5000000) 


rm(win, new_names)

# SIF data 


gps <- extract_gps('Uwe-3')[-c(1,2),]
flox <-
    read_delim('ALL_INDEX_FLOX_2019-06-28_09_21_26_locs.csv', delim = ',') %>%
    mutate(Date = substr(`datetime [UTC]`, 1,8),
           Time = as.POSIXct(`datetime [UTC]`, format = '%d-%m-%y %R')) %>%
    bind_cols(gps) %>%
    select(doy.dayfract,Time, Date, Location, Transect, contains("SIF"), x, y) %>%
    filter(Location != 0,
           x > 5) %>%
    group_by(Transect, Location) %>%
    summarise(doy.dayfract = mean(doy.dayfract),
              Date = mean(Time),
              SIF_A = mean(`SIF_A_ifld [mW m-2nm-1sr-1]`),
              SIF_B = mean(`SIF_A_ifld [mW m-2nm-1sr-1]`),
              SIF_A_sd = sd(`SIF_A_ifld [mW m-2nm-1sr-1]`),
              SIF_B_sd = sd(`SIF_A_ifld [mW m-2nm-1sr-1]`),
              Lat = mean(y),
              Long = mean(x))

locs <- center.of.gravity(flox$Long, flox$Lat, as.factor(flox$Location))
locs <- latlong_to_UTM(locs$x, locs$y, 32)

rm(gps)   

plotlist[["True Color Map with locs"]] <- 
    plotlist[["True Color Map"]] + 
    geom_point(data = locs, aes(X,Y, shape = as.factor(ID)), stroke = 2) + scale_shape_manual(values = c(1:10))


# plot the diurnal course of the reference location
plotlist[['2019-06-27 SIF reference']] <-
    read_delim('ALL_INDEX_FLOX_2019-06-28_09_21_26_locs.csv', delim = ',') %>%
    mutate(Date = substr(`datetime [UTC]`, 1,8),
           Time = as.POSIXct(`datetime [UTC]`, format = '%d-%m-%y %R')) %>%
    select(doy.dayfract,Time, Date, Location, Transect, contains("SIF")) %>%
    filter(Location == 1,
           Time > '2019-06-27') %>%
    ggplot(aes(Time, `SIF_A_ifld [mW m-2nm-1sr-1]`)) + geom_point()

# plot the diurnal course of 10 locations for two days of measurements

plotlist[['2019-06-27 SIF ~ Transect']] <- 
    flox %>% 
    filter(Date > '2019-06-27') %>% 
    ggplot(aes(Transect-6, SIF_A))  + geom_point() + 
    geom_segment(data = flox %>% 
                     filter(Date > '2019-06-27', Transect <12) %>% 
                     select(Transect, SIF_A) %>% 
                     bind_cols(flox %>% filter(Date > '2019-06-27', Transect >7) %>% select(Transect, SIF_A)),
                 aes(xend = Transect1 - 6,
                     yend = SIF_A1))

plotlist[['2019-06-27 SIF ~ Time']] <-
    flox %>%
    filter(Date > '2019-06-27') %>%
    ggplot(aes(Date, SIF_A, shape = as.factor(Location))) + 
        geom_point(size = 2, stroke = 1.5) + 
        scale_shape_manual(values = c(1:10)) +
        geom_errorbar(aes(ymin =SIF_A - SIF_A_sd, ymax = SIF_A + SIF_A_sd))


# reshape the original table to spread the locations as separate variables
flox2 <- 
    flox %>% 
    group_by(Transect) %>% 
    summarise(Date = mean(Date)) %>%
    select(Date) %>%
    bind_cols(flox %>% 
                  select(Location, SIF_A) %>% 
                  spread(Location, SIF_A)) %>% 
    filter(Date > '2019-06-27') 
    

# Calculate delta SIF for each possible combination
deltasif <- list()
cols <- c(3:12)

for(col in cols){
    deltasif[[col]] <- combn(pull(flox2, col), 2, diff)
}
rm(col, cols)
deltasif <- do.call(rbind,deltasif)

# bind with LAI 
deltasif <- as.data.frame(cbind((lai %>% filter(Day == '2019-06-27', Rep == 2) %>% select(LAI_mean))$LAI_mean, 
                                deltasif))

names(deltasif) <- c("LAI", c(1:15))



# Plot the delta sif against the lai 
plotlist[['2019-06-27 deltaSIF ~ LAI']] <- 
    deltasif %>% gather(Delta, deltaSIF, -LAI) %>% 
    ggplot(aes(LAI, deltaSIF, shape = Delta)) + 
    geom_point() + 
    scale_shape_manual(values=c(1:15))

# The above plot is too crowded
# calculate delta sif only between each consecutive transect 
flox3 <- 
    flox2 %>% 
    gather(Location, SIF_A, -c(Date, Transect)) %>% 
    select(Transect, SIF_A, Location) %>% 
    spread(Transect, SIF_A)

deltasif2 <- list()
for(i in c(3:7)){
    deltasif2[[i-2]] <- flox3[,i] - flox3[,i-1]
}
rm(i)
deltasif2 <- do.call(cbind, deltasif2)

# plot again with fewer deltas calculated 
plotlist[['2019-06-27 deltaSIF ~ LAI 2']] <- 
    deltasif2 %>% bind_cols(lai %>% 
                                filter(Day == '2019-06-27', Rep == 2) %>% 
                                select(LAI_mean)) %>% 
    gather(Delta, deltaSIF, -c(LAI_mean, Day, Rep)) %>% 
    ggplot(aes(LAI_mean, deltaSIF, shape = Delta)) + 
    geom_point() + geom_line()

plotlist[['2019-06-27 SIF ~ LAI ']] <- 
    flox3 %>%
        mutate(Location = as.numeric(Location)) %>% 
        left_join(lai, by = "Location") %>% 
        filter(Rep == 1,
               Day == '2019-06-27') %>% 
        gather(Window, SIF, -c(Location, Day, Rep, LAI_mean, LAI_sd)) %>%
        ggplot(aes(LAI_mean, SIF,  color = as.factor(Window))) + 
             # geom_line() + 
        geom_point(aes(shape = as.factor(Location)), stroke = stroke, size = size) + geom_smooth(method = 'lm', se = FALSE) +
             scale_shape_manual(values = 1:10) 
