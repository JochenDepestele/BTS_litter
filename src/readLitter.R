#################################################################################
## Helper functions for handling DATRAS litter exchange data
## Author: Casper W. Berg, DTU Aqua 
## Modifications for BTS analysis: David Vanavermaete and Jochen Depestele (ILVO)
#################################################################################

TVS2TVL <- function(){
    ## Swept area (m^2) for 30 min hauls
    TVS_TVL = c("TVLobs"=82798,"TVLobspred"=87163,"TVSobs"=68764,"TVSobspred"=68184)
    TVS_TVL
}

readlitter <- function (file = "IBTS.csv", na.strings = c("-9", "-9.0", "-9.00", "-9.0000"), type="Weight"){
  
  library(dplyr)

    getExt<-function(x) { l=nchar(x); substr(x,l-3,l) }
    if(getExt(file)==".zip"){
        tempdir <- tempdir()
        file <- unzip(file, exdir = tempdir)[1]
    }
    
    ## Locating lines with headers
    lines <- readLines(file)
    i <- grep("RecordType", lines)
    headers <- lines[i]
    skip <- i - 1
    nrow <- c(diff(i), 0) - 1
    cat("Reading data files\n")
    d <- lapply(1:length(i), function(i) {
        emptyDF <- function(nm) {
            structure(rep(list(logical(0)), length(nm)), .Names = nm, 
                class = "data.frame")
        }
        if (nrow[i] == 0) {
            return(emptyDF(strsplit(headers[i], ",")[[1]]))
        }
        
        ans <- read.csv(file, nrow = nrow[i], 
                                          skip = skip[i], na.strings = na.strings, stringsAsFactors = TRUE)
        
        ans$StNo <- as.character(ans$StNo)
        ans <- renameDATRAS(ans)
        ans
    })
    names(d) <- sapply(d, function(x) {
        as.character(x$RecordType[1])
    })

    if (sum(sapply(d, nrow)) + length(i) != length(lines)) 
        stop("csv file appears to be corrupt.")
    if (is.null(d[[1]]$StatRec)) 
        d[[1]]$StatRec <- d[[1]]$AreaCode
    haul.id <- quote(factor(paste(Year, Quarter, Country, Ship, 
                                  Gear, StNo, HaulNo, sep = ":")))

    d[[1]]$haul.id = eval(haul.id,d[[1]])
    d[[2]]$haul.id = eval(haul.id,d[[2]])
    d[[2]]$haul.id = factor(d[[2]]$haul.id, levels=levels(d[[1]]$haul.id))
    
    LTnotHH = setdiff(levels(d[[2]]$haul.id),levels(d[[1]]$haul.id))
    HHnotLT = setdiff(levels(d[[1]]$haul.id),levels(d[[2]]$haul.id))
    if( length(LTnotHH) > 0 ){
        print(paste("Dropping",length(LTnotHH),"LT records with no HH info:"))
        print(LTnotHH)
        d[[2]] = subset(d[[2]],!haul.id %in% LTnotHH)
    }

    if( length(HHnotLT) > 0 ){
        print(paste("Dropping",length(HHnotLT),"HH records with no LT info:"))
        print(LTnotHH)
        d[[1]] = subset(d[[1]],!haul.id %in% HHnotLT)
    }

    
    ## Zero litter hauls are recorded as RECO-LT and PARAM = LT-TOT...don't exclude!
    badhauls = unique(d$LT$haul.id[ d$LT$LTREF=="RECO-LT" & d$LT$PARAM!="LT-TOT" ] )
    if(length(badhauls)>0){
        cat("Dropping ", length(badhauls), " hauls of type RECO-LT with PARAM!=LT-TOT\n")
        d[[2]] = subset(d[[2]], !haul.id %in% badhauls)
        d[[1]] = subset(d[[1]], !haul.id %in% badhauls)
    }
    
        
    if(type=="Numbers"){
        badhauls = unique(d$LT$haul.id[ is.na(d$LT$LT_Items) ] )
        badhauls <- c(badhauls, # Three hauls (1 haul in 2017, 2018 and 2019) have unrealistic high number of litter items compared to their weight in grammes. They are excluded
                      d$LT$haul.id[ d$LT$LT_Weight<1000 & d$LT$UnitWgt=="g/haul" & d$LT$LT_Items>100 ])
        if(length(badhauls)>0){ 
            cat("Dropping ", length(badhauls), " hauls with only weight info\n")
            d[[2]] = subset(d[[2]], !haul.id %in% badhauls)
            d[[1]] = subset(d[[1]], !haul.id %in% badhauls)
        }
    }
    ## zero hauls: weight unit does not matter
    cat("Weight units:\n")
    print(levels(d$LT$UnitWgt))
    d$LT$UnitWgt <- factor(d$LT$UnitWgt,levels=c("g/haul","kg/haul"))
    
    zerohauls = d$LT$LT_Items==0 | d$LT$LT_Weight==0
    d$LT$UnitWgt[ zerohauls ] = "kg/haul"
        
    if(type=="Weight"){
      d$LT$LT_Weight[ !is.na(d$LT$LT_Weight) & d$LT$LT_Weight>100 & d$LT$UnitWgt=="kg/haul" & d$LT$Ship=="11BE" & d$LT$LT_Items==1] # 4 hauls in DATRAS with incorrect weight (Heleen Raat - Belgian BTS coordination, pers comm 14 Feb 2023)
      d$LT$UnitWgt[ !is.na(d$LT$LT_Weight) & d$LT$LT_Weight>100 & d$LT$UnitWgt=="kg/haul" & d$LT$Ship=="11BE" & d$LT$LT_Items==1] <- "g/haul"
      d$LT$LT_Weight[ !is.na(d$LT$LT_Weight) & d$LT$LT_Weight>100 & d$LT$UnitWgt=="kg/haul" & d$LT$LT_Items==1]
      badhauls = unique(d$LT$haul.id[ is.na(d$LT$UnitWgt) | !d$LT$UnitWgt%in%c("g/haul","kg/haul")])
      badhauls = c(badhauls,
                   d$LT$haul.id[ !is.na(d$LT$LT_Weight) & d$LT$LT_Weight>100 & d$LT$UnitWgt=="kg/haul" & d$LT$LT_Items==1])
        if(length(badhauls)>0){ 
            cat("Dropping ", length(badhauls), " hauls with no or bad weight unit info\n")
            d[[2]] = subset(d[[2]], !haul.id %in% badhauls)
            d[[1]] = subset(d[[1]], !haul.id %in% badhauls)   
        }
        d$LT$LT_Weight[ d$LT$UnitWgt=="g/haul" ] = d$LT$LT_Weight[ d$LT$UnitWgt=="g/haul" ] / 1000
    }
    
    if(nlevels(d[[1]]$haul.id) != nlevels(d[[2]]$haul.id)){
        print("Missing haul info:")
        cat("HH: " ,nlevels(d$HH$haul.id)," ,LT: ",nlevels(d$LT$haul.id),"\n","Setdiff:\n")
        print( setdiff(levels(d[[1]]$haul.id),levels(d[[2]]$haul.id)))
        stop()
    }
    
    conv = read.csv("./data/litter_conversion.csv")
    lookup <- function(x,table){
        x <- factor(x)
        levels(x) <- table[levels(x)]
        if(is.numeric(table)) as.numeric(as.character(x)) else x
    }

    Code2Type <- structure(as.character(conv$Type), names=conv$C.TS)
    Code2TypeRev <- structure(as.character(conv$Type), names=conv$C.TS.REV)
    d[[2]] = transform(d[[2]],Type= factor(
                      ifelse(as.character(LTREF)=="C-TS",
                             as.character( lookup(PARAM,Code2Type) ),
                             as.character( lookup(PARAM,Code2TypeRev) )
                             )))

    d[[2]]$Type <- factor(d[[2]]$Type,levels=c("Glass","Metal","Natural","Other","Plastic","Rubber"))
                      
    if(type=="Weight"){
        agg = xtabs(LT_Weight ~ haul.id + Type,d[[2]])[ d[[1]]$haul.id, ]
        
    } else if(type=="Numbers"){
        agg = xtabs(LT_Items ~ haul.id + Type,d[[2]])[ d[[1]]$haul.id, ]
    }
    dn = dimnames(agg)[2][[1]]
    for(i in 1:length(dn))
        d[[1]][,dn[i]] = agg[,i]

    ## Extra categories:
    ## 1. Singe use plastic (SUP)
    Code2Type <- structure(as.character(conv$SUP), names=conv$C.TS)
    Code2TypeRev <- structure(as.character(conv$SUP), names=conv$C.TS.REV)
    d[[2]] = transform(d[[2]],SUP= factor(
                                  ifelse(as.character(LTREF)=="C-TS",
                                         as.character( lookup(PARAM,Code2Type) ),
                                         as.character( lookup(PARAM,Code2TypeRev) )
                                         )))
    if(type=="Weight"){
        agg = xtabs(LT_Weight ~ haul.id + SUP,d[[2]])[ d[[1]]$haul.id, ]
        
    } else if(type=="Numbers"){
        agg = xtabs(LT_Items ~ haul.id + SUP,d[[2]])[ d[[1]]$haul.id, ]
    }
    d[[1]]$SUP = agg[ ,"Yes"]  
         
    ## 2. Fishing related
    Code2Type <- structure(as.character(conv$Fishing.related), names=conv$C.TS)
    Code2TypeRev <- structure(as.character(conv$Fishing.related), names=conv$C.TS.REV)
    d[[2]] = transform(d[[2]],Fishing.related = factor(
                                  ifelse(as.character(LTREF)=="C-TS",
                                         as.character( lookup(PARAM,Code2Type) ),
                                         as.character( lookup(PARAM,Code2TypeRev) )
                                         )))
    if(type=="Weight"){
        agg = xtabs(LT_Weight ~ haul.id + Fishing.related,d[[2]])[ d[[1]]$haul.id, ]
        
    } else if(type=="Numbers"){
        agg = xtabs(LT_Items ~ haul.id + Fishing.related,d[[2]])[ d[[1]]$haul.id, ]
    }
    d[[1]]$Fishing.related = agg[,"Yes"]

    ## Extra variables (for surveyIndex package)
    d[[1]]$lon = (d[[1]]$ShootLong + d[[1]]$HaulLong)/2
    d[[1]]$lat = (d[[1]]$ShootLat + d[[1]]$HaulLat)/2
    
    d[[1]]$ctime = d[[1]]$Year + (d[[1]]$Month-1)/12 + (d[[1]]$Day-1)/365

    d[[1]]$timeOfYear = (d[[1]]$Month-1)/12 + (d[[1]]$Day-1)/365
    
    d[[1]]$Year = factor(d[[1]]$Year)

    d[[1]]$TimeShotHour = 12

    # in case we'd like to remove outliers:    
    # ## check Litter data on inconsistency
    # temp<-d[[2]]%>%group_by(haul.id)%>%summarize(LT_Items)
    # ll<-mean(temp$LT_Items, na.rm=T)-3*sd(temp$LT_Items, na.rm=T)
    # ul<-mean(temp$LT_Items, na.rm=T)+3*sd(temp$LT_Items, na.rm=T)
    # 
    # outliers<-temp$haul.id[which(temp$LT_Items>ul)]
    # d[[2]]%>%filter(!(haul.id%in%outliers)) -> d[[2]]
    # d[[1]]%>%filter(!(haul.id%in%outliers))  -> d[[1]]
    
    #incorrect reported
    excl<-d[[2]] %>% filter((LT_Items ==0 & PARAM!="LT-TOT")) %>% pull(haul.id)
    

    d[[2]] %>% filter(!(haul.id%in%excl)) -> d[[2]]
    d[[1]] %>% filter(!(haul.id%in%excl)) -> d[[1]]

    
    
    ## check distance
    r_earth<-6378.135 #km
    L_s<-d[[1]]$ShootLat
    l_s<-d[[1]]$ShootLong
    L_h<-d[[1]]$HaulLat
    l_h<-d[[1]]$HaulLong
    
    
    D= 1000* r_earth * acos ( cos( (90-L_s) * (pi/180)) * cos((90-L_h)* (pi/180))  + 
         sin((90-L_s)*(pi/180))*sin((90-L_h)*(pi/180))*cos((l_s-l_h)* (pi/180)) )
    
    d[[1]]$Distance_Calculated<-D
    D_used<-d[[1]]$Distance
    
    #if Distance is NA used calculated distance
    D_used[which(is.na(D_used))]<- D[which(is.na(D_used))]
    
    D_used[which( D*0.8>=d[[1]]$Distance | D*1.2<= d[[1]]$Distance)]<- D[which( D*0.8>=d[[1]]$Distance | D*1.2<= d[[1]]$Distance)]
    d[[1]]$Distance_Used<-D_used

    TVS_TVL = TVS2TVL()
    
    #isolate the beam length
    beam_width<-as.numeric(gsub(".*?([0-9]+).*", "\\1", as.character(d[[1]]$Gear)))
    
    #Surface = distance of haul times beam length
    d[[1]]$EFFORT = d[[1]]$Distance_Used * beam_width
    
    #d[[1]]$EFFORT = TVS_TVL[2]/30 * d[[1]]$HaulDur
    #d[[1]]$EFFORT[ d[[1]]$Gear == "TVS" ] = TVS_TVL[4]/30 * d[[1]]$HaulDur[ d[[1]]$Gear == "TVS" ]
    
    # estimate EFFORT in m?? !!!
    
    d[[1]]
}



mybubblePlot<-function (d, response = "HaulWgt", scale = NULL, col.zero = "red", 
    pch.zero = "+", rim = FALSE, ...) 
{
    d$resp.var <- d[,response]
    if (is.null(scale)) 
        scale = 4/max(sqrt(d$resp.var), na.rm = TRUE)
    plot(d$lon, d$lat, type = "n", xlab = "Longitude", ylab = "Latitude",...)
    map("worldHires", fill = TRUE, plot = TRUE, add = TRUE, col = grey(0.5))
    points(d$lon, d$lat, pch = 16, cex = scale * sqrt(d$resp.var))
    if (rim) 
        points(d$lon, d$lat, cex = scale * sqrt(d$resp.var), 
            pch = 1, lwd = 0.5, col = "blue")
    zero = subset(d, resp.var == 0)
    points(zero$lon, zero$lat, pch = pch.zero, col = col.zero)
}
