correrEjemplo=function(){
  #La función ejecuta el caso de prueba Brompton de DynaTopmodel
  #Como se encontró un error en la función run.dtm() del paquete
  #se usó la función del script que se renombró run.dtm2()
  
  library(dynatopmodel)
  library(zoo)
  library(xts)
  
 
  data(brompton)
  sp::plot(brompton$dem)
  sp::plot(brompton$drn)
  
  chans <- build_chans(dem=brompton$dem, drn=brompton$drn, chan.width=2)  #buff=5, 
  sp::plot(chans)
  
  sp::plot(brompton$flowdists)
  # discretisation by reverse distance from nearest channel. The raster brompton$flowdists 
  # gives the D8 flow pathway distance for every area in the catchment  
  layers <- addLayer(brompton$dem, 2000-brompton$flowdists) 
  sp::plot(layers)
  
  disc <- discretise(layers, cuts=c(flowdists=5), chans=chans, area.thresh=3/100)
  write.table(disc$groups, sep="\t", row.names=FALSE)
  
  #Network routing table 
  routing <- build_routing_table(brompton$dem, chans)
  
  # Here we apply the same parameter values to all groups. Suggest applying smaller m and td values to 
  # the closest areas to simulate a fast response due to the artificial drainage. 
  # It would also be possible to supply a custom transmissivity profile that has 
  # a discontinuity at the depth of the drains 
  groups <- disc$groups 
  groups$m <- 0.011 
  groups$td <- 42 
  # a very high transmissivity prevents saturation flow as there appears be little 
  groups$ln_t0 <- 18 
  groups$srz_max <- 0.1 
  # initial root zone storage 
  groups$srz0 <- 0.87 
  # quite slow channel flow, which might be expected with the shallow and reedy 
  # reaches in this catchment 
  groups$vchan <- 750
  
  
  # Observations at a 15 minute time step 
  dt <- 0.25 
  obs <- list(rain=brompton$rain, pe=brompton$pe, qobs=brompton$qobs) 
  obs <- aggregate_obs(obs, dt=dt)
  
  
  # parameters for graphics output 
  par <- disp.par(int.time=24)
  # Note max.q in mm/hr 
  par$max.q <- 1000*max(obs$qobs, na.rm=TRUE) 
  sim.start <- "2012-09-23"
  sim.end <- "2012-10-01"
  
  # Ensure output goes to a new window 
  options("device"="X11") 
  # take initial discharge from the observations 
  qt0 <- as.numeric(obs$qobs[sim.start][1])
  

  # Run the model across the September 2012 storm event using 2 inner time steps and a 15 minute interval 
  setwd("C:/Users/mariyeg/Documents/Dynatopmodel/lib")
  source('list_util.R')
  source('calib_util.R')
  source('disp_util.R')
  source('evap.R')
  source('hydro_util.R')
  source('init.R')
  source('maths_util.R')
  source('plot_util.R')
  source('rast_util.R')
  source('read_obs.R')
  source('xts_util.R')
  
  setwd("C:/Users/mariyeg/Documents/Dynatopmodel/model")
  source('balance.R')
  source('debug.R')
  source('defs.R')
  source('distribute.R')
  source('dtm_main.R')
  source('initial.R')
  
  #debug(init.input) 
  
  source('kinematic_desolve.R')
  source('model_results.R')
  source('rain.R')
  source('root_zone.R')
  source('route_ovf.R')
  source('routing.R')
  source('unsat_zone.R')
  source('update_subsurface.R')
  
  storm.run <- run.dtm2(groups=groups, weights=disc$weights, rain=obs$rain, 
                       pe=obs$pe, qobs=obs$qobs, qt0=qt0, sim.start=sim.start,
                       sim.end=sim.end, routing=routing, disp.par=par, ntt=2)
  # show run statistics 
  cat("NSE=", NSE(storm.run$qsim, storm.run$qobs)) 
  cat("Time at peak =", format(time_at_peak(storm.run$qsim)))
  
}

correrEjemplo()