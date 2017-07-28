leyendoDem=function(){
  #Esta función leer el DEM raster y define la extensión espacial
  #También reclasifica para eliminar valores no data
  
  library(lattice)
  library(sp)
  library(raster)
  
  midir<-getwd()
  nuevodir<-paste(midir, "/mapas", sep = "")
  setwd(nuevodir)
  DEM<-raster(x="Valle_derecho_dem.mpr")
  setwd(midir)
  
  # Definiendo xMin,xMax,ymin,yMax
  #Extensión para resolución 25m
  #rasExt <- extent(291579.5, 304454.5, 962403, 988178)
  #Extensión para resolución 5m
  rasExt <- extent(291579.5, 294154.5, 962398, 967553)
  
  DEM@extent <- rasExt
  #Necesario reclasificar para sustituir por NA lo que no es parte del mapa
  DEMR <- reclassify(DEM, c(-Inf,3600,NA))
  print(DEMR)
  print(summary(DEMR))
  plot(DEMR,main="Modelo de Elevación Digital")
  
  return(DEMR)
  
}