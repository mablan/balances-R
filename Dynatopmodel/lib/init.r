require(utils)
require(tools)

# mirroring the HEC file structure
#scratch.dir <- "/scratch/hpc/30/metcalfp"
#source("~/source/set_paths.r") #MARILENA


# default location
# if(!file.exists(Sys.getenv("SRC_DIR")))
# {
#   Sys.setenv("SRC_DIR"=path.expand("~/source"))
#   cat("Setting source directory to ", Sys.getenv("SRC_DIR"), "\n")
# }

#source("~/source/load_libs.r") #MARILENA

# LoadPackages(pnames=c("fields",
# #"rgl", "dismo", "plotKML",
# "akima", "intervals", "raster", "spam",
# "sp", "deSolve",
# "lattice", "rgdal",
# "maptools", "igraph", "rgeos"))

require("fields")
require("akima") 
require("intervals")
require("raster") 
require("spam")
require("sp")
require("deSolve")
require("lattice")
require("rgdal")
require("maptools")
require("igraph")
require("rgeos")

# LoadPackages(pnames=c("xts", "raster", "lattice", "rgdal", "deSolve",
#                       "stringr", "lubridate", "igraph", "rgeos", "topmodel", "shape"), attach=T)

require("xts")
require("raster")
require("lattice") 
require("rgdal")
require("deSolve")
require("stringr")
require("lubridate")
require("igraph")
require("rgeos")
require("topmodel")
require("shape")


#options("stringsAsFactors"=F)

# force default device to windows (stops RStudio from only only its plot window to be used)
options(device = "X11")

# useful CRS
osgb <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +towgs84=446.448,-125.157,542.06,0.15,0.247,0.842,-20.489 +units=m +no_defs"
# # spherical mercator
# get.epsg(latlong) = "+init=epsg:4326"
lonlat = "+init=epsg:4326"
# ukgrid = "+init=epsg:27700"
google = "+init=epsg:3857"