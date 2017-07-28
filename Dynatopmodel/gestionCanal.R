crearCanales=function() {
  #Esta función ayudo a correr TOPMODEL 
  #especificamente a construir la matriz delay
  
  source('gestionDem.R')
  DEMR<-leyendoDem()  
  library(dynatopmodel)
  
  #Se determinan las áreas de contribución pendiente arriba (Es un wrapper de la función implementada en TOPMODEL, topidx)
  #Opcionalmente calcula el indice de humedad topografica
  areaindice<- upslope.area(DEMR, atb=TRUE, deg = 0.1, fill.sinks = TRUE)
  #areaindice (Raster Stack, layers: list of 2)
  sp::plot(areaindice, main=c("Upslope area (log(m^2/m))", "TWI log(m^2/m)"))
  
  #Construye raster para la localización del canal, a  partir del indice de humedad topográfica
  canales<-build_chans(atb = areaindice$atb, chan.width = 1, atb.thresh = 0.5)
  #canales raster de 2 bandas (Raster Stack, layers: list of 2)
  #la primera celdas con canales y la segunda proporciones
  sp::plot(canales[[1]], maxpixels=400000, main="Canales")
  
  #Genera una tabla de ancho de red para la cuenca
  #Cuando se pasa a la rutina run.dtm se usa para enrrutar
  #flujo de canal al aforo durante la corrina de DynamicTopmodel
  enrrutamiento <- build_routing_table(DEMR, chans=canales, breaks=4)
  #A two-column data.frame. Its first column is the average flow distance to the outlet, in m, the second
  #the proportions of the catchment channel network within each distance category
  print(enrrutamiento)
  return(canales)
}

crearCanales()