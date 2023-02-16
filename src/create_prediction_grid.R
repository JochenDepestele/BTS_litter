#################################################################################
## # Get prediction grid with bathymetry data using EMODnet-bathymetry
## Author: Jochen Depestele (ILVO)
#################################################################################
library(dplyr)
library(sf)
library(raster)
library(sp)
library(raster)
library(mapplots)
# library(sfdSAR)
library(vmstools)
library(data.table)

OneDrivedir <- file.path(Sys.getenv("OneDrive")) # This script should find the location of your OneDrive directory
datadir <- paste0(OneDrivedir,"/000_data")
utildir <- paste0(OneDrivedir,"/000_utilities")

source(paste0(file.path(utildir),"/coords_to_csquare_VMStools.R"))


ICESareas <- sf::st_read(dsn = "./shapefiles/ICES/ICES_areas.shp")
sf_squares <- st_as_sf(st_make_grid(ICESareas,cellsize = c(0.05,0.05), # use st_as_sf to convert sfc_POLYGON to sf data.frame
                                    offset = st_bbox(ICESareas)[c("xmin", "ymin")],
                                    square = TRUE,flat_topped = FALSE))

ROI <- extent(st_bbox(st_transform(ICESareas,crs=4326)))
# Get EMODnet bathymetry data from the nc files ####
d4 <- paste0(file.path(datadir),"/EMODnet-Bathymetry","/D4_2020.nc") # Downloaded from https://portal.emodnet-bathymetry.eu/# (downloaded 10 October 2022)
d4 <- brick(d4,varname="elevation")
d5 <- paste0(file.path(datadir),"/EMODnet-Bathymetry","/D5_2020.nc") # Downloaded from https://portal.emodnet-bathymetry.eu/# (downloaded 10 October 2022)
d5 <- brick(d5,varname="elevation")
e3 <- paste0(file.path(datadir),"/EMODnet-Bathymetry","/E3_2020.nc") # Downloaded from https://portal.emodnet-bathymetry.eu/# (downloaded 10 October 2022)
e3 <- brick(e3,varname="elevation")
e4 <- paste0(file.path(datadir),"/EMODnet-Bathymetry","/E4_2020.nc") # Downloaded from https://portal.emodnet-bathymetry.eu/# (downloaded 10 October 2022)
e4 <- brick(e4,varname="elevation")
e5 <- paste0(file.path(datadir),"/EMODnet-Bathymetry","/E5_2020.nc") # Downloaded from https://portal.emodnet-bathymetry.eu/# (downloaded 10 October 2022)
e5 <- brick(e5,varname="elevation")
raster::crs(d4) <- "EPSG:4326"
raster::crs(d5) <- "EPSG:4326"
raster::crs(e3) <- "EPSG:4326"
raster::crs(e4) <- "EPSG:4326"
raster::crs(e5) <- "EPSG:4326"
d4_crop <- crop(d4,ROI) # Note: if mask= F, the crop will be by extent (box) 
d5_crop <- crop(d5,ROI)
e3_crop <- crop(e3,ROI)
e4_crop <- crop(e4,ROI)
e5_crop <- crop(e5,ROI)
bathy <- merge(d4_crop, d5_crop, e3_crop, e4_crop, e5_crop)
plot(bathy)
# plot(ICESareas,add=T)
ICESareas <- st_transform(ICESareas,crs=4326)
sp_ICESareas <- as(ICESareas,'Spatial')
sp_ICESareas$ICES_SUB <- as.factor(sp_ICESareas$ICES_SUB)
r_ICESareas <- raster::rasterize(sp_ICESareas, bathy, 
                                 field=sp_ICESareas$ICES_SUB, method="ngb")
plot(r_ICESareas)
bathy <- mask(bathy,r_ICESareas)
plot(bathy)
# sf_squares <- st_transform(sf_squares,crs=4326)

sf_squares <- st_transform(sf_squares, crs=27700) # st_centroid requires planar coordinates to calculate distances
sf_squares <- sf_squares %>%
  dplyr::mutate(centroid = sf::st_centroid(.))
sf_squares$centroids <- st_transform(sf_squares, 27700) %>% 
  st_centroid() %>% 
  st_transform(., 4326) %>%
  st_geometry()
sf_squares <- sf_squares %>% 
  dplyr::mutate(lon_csq = sf::st_coordinates(centroids)[,1], # Add midpoint longitude of centroids
                lat_csq = sf::st_coordinates(centroids)[,2], # Add midpoint latitude of centroids
                c_square = CSquare(lon_csq,lat_csq,degrees=0.05),
                stat_rec = ices.rect2(lon_csq, lat_csq))

lonlat <- as.data.frame(sf_squares %>% dplyr::select(lon_csq, lat_csq) %>% st_drop_geometry())
sf_squares$depth_m <- abs(raster::extract(bathy, lonlat))

# Joining the ICESareas dataframe with the sf_squares retains many gridcells on land
ICESareas_dt <- setDT(get(load(paste0(file.path(datadir),"/Areas/ICES/rICES/ICESAreaRects.rdata"))))
selcols <- c("AreaName","StatRect","subdivision")
ICESareas_dt <- unique(ICESareas_dt[,..selcols])
sf_squares <- sf_squares %>%
  left_join(ICESareas_dt,by=c('stat_rec'='StatRect'))

bgrid <- sf_squares %>%
  dplyr::select(lon_csq,lat_csq,c_square,stat_rec,depth_m,subdivision) %>%
  st_drop_geometry()

bgrid <- bgrid[!is.na(bgrid$subdivision),]

names(bgrid) <- c("lon","lat","c_square","stat_rec","Depth","ICES_SUB")

unique(bgrid[which(bgrid$ICES_SUB %in% c("IVa","IVb","IVc","VIId","VIIa","VIIf","VIIg")),]$ICES_SUB)
nrow(bgrid[which(bgrid$ICES_SUB %in% c("IVa","IVb","IVc","VIId","VIIa","VIIf","VIIg")),])

write.csv2(bgrid[which(bgrid$ICES_SUB %in% c("IVa","IVb","IVc","VIId","VIIa","VIIf","VIIg")),],file=paste0("./data/bgrid.csv"),row.names = F)


