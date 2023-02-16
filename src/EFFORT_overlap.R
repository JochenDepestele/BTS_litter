library(DATRAS)
library(mgcv)
library(nlme)
library(surveyIndex)
library(data.table)
library(sf)


# Load the prediction grid
bgrid <- read.csv2(paste0("C:/Users/jdepestele/OneDrive - ILVO/gitr/SEAwise/SEAwise_T4-5/BTS_litter/data/bgrid.csv"))
bgrid <- bgrid[,c("lon","lat","Depth","ICES_SUB","c_square")]
bgrid = subset(bgrid, ICES_SUB %in% c("VIIa","VIIf","VIIg"))
nrow(bgrid) # 9200 grid cells

# Load fishing effort
OneDrivedir <- file.path(Sys.getenv("OneDrive")) # This script should find the location of your OneDrive directory
datadir <- paste0(OneDrivedir,"/000_data")
sf_fishing <- readRDS(paste0(datadir,"/OSPAR_SAR/fishing_Rijnsdorp2020metiers_20092020.RDS"))
fishing <- sf_fishing %>% st_drop_geometry() %>% setDT()

# Load models and their predictions
models <- readRDS("C:/Users/jdepestele/OneDrive - ILVO/gitr/SEAwise/SEAwise_T4-5/BTS_litter/output/NWW/models.RDS")
models2 <- readRDS("C:/Users/jdepestele/OneDrive - ILVO/gitr/SEAwise/SEAwise_T4-5/BTS_litter/output/NWW/models2.RDS")
models3 <- readRDS("C:/Users/jdepestele/OneDrive - ILVO/gitr/SEAwise/SEAwise_T4-5/BTS_litter/output/NWW/models3.RDS")

y = which(as.numeric(as.character(names(models[["Fishing.related"]]$gPreds2))) == year)
y = models[["Fishing.related"]]$yearNum
surveyIndex:::concTransform(log(models[["Fishing.related"]]$gPreds2))

nrow(models[["Fishing.related"]]$gPreds[[1]])
class(models[["Fishing.related"]]$gPreds2)
length(models[["Fishing.related"]]$gPreds2[[1]])
names(models[["Fishing.related"]]$gPreds2[[1]])

nrow(models[["Fishing.related"]]$gPreds2[[1]][["2015"]]) # Predictions for 2015

