library(rgdal)
library(tmap)
library(maptools)
library(tmap)
library(spatstat)
library(gstat) # Use gstat's idw routine
library(sp)    # Used for the spsample function
library(raster)
library(rgeos)
setwd("D:/ancestry/India Shape/")
d<-read.csv("lat_long",sep="\t", header=TRUE)


###
UTM32n <- CRS("+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
# World Geographic System 1984 (lat/long) - mapping
WGS84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 
p <-  SpatialPointsDataFrame(coords = d[,c("Longitude", "Latitude")], 
                                  data = d, 
                                  proj4string = UTM32n)

# Load India boudary map
shp<-readOGR("D:/ancestry/India Shape/IGISMAP/Indian_States.shp")

# Replace point boundary extent with that of India
p@bbox <- shp@bbox

tm_shape(shp) + tm_polygons() +
  tm_shape(p) +
  tm_dots(col="G", palette = "YlOrRd",
          title="rs1799971", size=0.7) +
  tm_text("G", just="left", xmod=.5, size = 0.7) +
  tm_legend(legend.outside=TRUE)
####op warning: Warning messages:
#1: The argument auto.palette.mapping is deprecated. Please use midpoint for numeric data and stretch.palette for categorical data to control the palette mapping. 
#2: Currect projection of shape shp unknown. Long-lat (WGS84) is assumed. 

th  <-  as(dirichlet(as.ppp(p)), "SpatialPolygons")

# The dirichlet function does not carry over projection information
# requiring that this information be added manually
proj4string(th) <- proj4string(p)

# The tessellated surface does not store attribute information
# from the point data layer. We'll use the over() function (from the sp
# package) to join the point attributes to the tesselated surface via
# a spatial join. The over() function creates a dataframe that will need to
# be added to the `th` object thus creating a SpatialPolygonsDataFrame object
th.z     <- over(th, p, fn=mean)
th.spdf  <-  SpatialPolygonsDataFrame(th, th.z)

# Finally, we'll clip the tessellated  surface to the Texas boundaries
set_RGEOS_CheckValidity(2L)
th.clp   <- intersect(shp,th.spdf)
### x[subsx, ] is invalid
#Error in RGEOSBinTopoFunc(spgeom1, spgeom2, byid, id, drop_lower_td, unaryUnion_if_byid_false,  : 
#                            TopologyException: Input geom 0 is invalid: Ring Self-intersection at or near point 79.775711259999994 11.83757306 at 79.775711259999994 11.83757306
#                          In addition: Warning messages:
#                            1: In RGEOSUnaryPredFunc(spgeom, byid, "rgeos_isvalid") :
#                            Ring Self-intersection at or near point 79.775711259999994 11.83757306
#                          2: In rgeos::gIntersection(x[subsx, ], y[subsy, ], byid = TRUE, drop_lower_td = TRUE) :
#                            Invalid objects found; consider using set_RGEOS_CheckValidity(2L)

# Map the data
tm_shape(th.clp) + 
  tm_polygons(col="G", palette="YlOrRd", auto.palette.mapping=FALSE,
              title="rs1799971") +
  tm_legend(legend.outside=TRUE)

grd              <- as.data.frame(spsample(p, "regular", n=50000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object

# Add P's projection information to the empty grid
proj4string(p) <- proj4string(p) # Temp fix until new proj env is adopted
proj4string(grd) <- proj4string(p)

# Interpolate the grid cells using a power value of 2 (idp=2.0)
P.idw <- gstat::idw(G ~ 1, p, newdata=grd, idp=2.0)

# Convert to raster object then clip to Texas
r       <- raster(P.idw)
r.m     <- mask(r, shp)

# Plot
tm_shape(r.m) + 
  tm_raster(n=10,palette = "YlOrRd", auto.palette.mapping = FALSE,
            title="rs1799971") + 
  tm_shape(p) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)+tm_legend(legend.outside=TRUE)+tm_text("POP", just="top", xmod=0.7, size = 0.6)
