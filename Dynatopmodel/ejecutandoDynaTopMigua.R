ejecutarDyna=function() {
  #La función intenta ejecutar Dynatopmodel para Miguaguó
  source('gestionDem.R')
  source('gestionCanal.R')
  
  DEMR<-leyendoDem() 
  canales<-crearCanales()
 
  
  #Pensando en discretizar se calcula la pendiente
  #terrain del paquete Raster, when neigbors=8, slope is computed according to Horn (1981).
  #The Horn algorithm may be best for rough surfaces
  tangente<-terrain(DEMR, opt = 'slope', unit='tangent')
  #pendiente en %
  pendiente<-tangente*100
  print(summary(pendiente))
  hist(pendiente)
  sp::plot(pendiente, main="Pendiente en %")
  
  #Más o menos como propone córdova
  #unidades<- reclassify(pendiente, c(0,2,1,2,5,2,5,15,3,15,46,4))
  #capas<- addLayer(DEMR, unidades)
  #sp::plot(capas, main=c("DEM", "UNIDADES"))
  capas<- addLayer(DEMR, pendiente)
  sp::plot(capas, main=c("DEM", "PENDIENTES"))
  
  discretiza <- discretise(capas, cuts=c(slope=5), chans=canales, area.thresh=3/100)
  write.table(discretiza$groups, sep=",", row.names=FALSE)
  #par(mfrow=c(1,1))
  sp::plot(discretiza$hru, main="HRU", col=c("red","blue","green"))
  hist(discretiza$hru, label=TRUE)
  
}

ejecutarDyna()