#source("kinematic_desolve.r")
#c.route.kinematic.euler <- cmpfun(route.kinematic.euler)
################################################################################
# run inner loop Solve kinematic equation for base flow at given base flows at
# previous steps and input at this one
# Input (all)
# w               : weighting matrix
# pe              : potential evap
# tm,             : current simulation time
# ntt,            : no. inner time steps
# dt,             : main time step
# dqds,           : gradient functions
# ------------------------------------------------------------------------------
# Input (each HSU)
# ------------------------------------------------------------------------------
# flows$qbf       : base (specific) flow  at previous time step
# flows$qin       : input (total) flow at previous time step
# flows$uz        : drainage recharge from unsaturated zone (m/hr)
# groups$area : plan area
#
# stores$sd       : specific storage deficits
# ichan           : channel identifiers
# ------------------------------------------------------------------------------
# Input
# ------------------------------------------------------------------------------
# Returns (each unit)
# ------------------------------------------------------------------------------
# flows$qbf       : specific base flow at time step, per grouping
# flows$qin       : total input (upslope) input, per grouping
# stores$sd       : updated specific storage deficits: unsat drainage and qin inputs, qbf out
# flows$qof       : updated specific overland flux per areal grouping
################################################################################
update.subsurface <- function (groups, flows, stores,
						              w,           # weighting matrix
						              A,           # complementary weighting matrix A= diag(1/a, N, N) %*% t(w) %*% diag(a, N, N) - identity.matrix(N)
                          pe=0,        # potential evap
                          tm,          # current simulation time
                          ntt,         # no. inner time steps
                          dt,          # main time step
                          dqds=NULL,   # gradient functions
                          ichan=1,
						              log=NULL)
{
  dtt <- dt/ntt

  # record of actual evapotranpiration and unsaturated drainage over each inner step
  ae.dtt <- matrix(0,ncol=nrow(groups), nrow=ntt)
  uz.inner <- ae.dtt 
  
  timei<-tm
  
  # if(tm == as.POSIXct("2015-11-20 09:30:00"))
  # {
  #   browser()
  # }
  # total input into channel
  Qriv <- rep(0, length(ichan))
  
  # area of the channel network
  achan <- groups[ichan,]$area
  
  flows$ex <- 0

	ex <- 0
	
	a <- groups$area 
	N <- nrow(w)
#	A <- 
	stores0 <- stores

  for(inner in 1:ntt)
  {
  	# save storages
  	sd0 <- stores$sd
  	
  	# apply rain input and evapotranspiration (*rates*) across the inner time step
    # note that rain and actual evap are maintained in flows
    updated <- root.zone(groups, flows, stores, pe,
                         dtt, timei, ichan)        
    
    # ae is removed from root zone only - note that original DynaTM has code to remove evap
    # at max potential rate from unsat zone
    flows <- updated$flows  # includes ae and rain
  	stores <- updated$stores  #  storage excess

    # route excess flow from root zone into unsaturated zone and  drainage into water table
    updated <- unsat.zone(groups, flows, stores, dtt, timei, ichan)
  	flows <- updated$flows
    stores <- updated$stores
    
    # Distribute baseflows downslope through areas using precalculated inter-group
    # solution of ODE system. Uses the Livermore solver lsosa
    flows <- route.kinematic.euler(groups, flows, stores, dtt,
    								ichan=ichan, w=w, time=timei, A=A,
                    dqds=dqds)
    
    # distribute fluxes to give new estimate for qin(t) given qbf(t) determined above
    # base flux transferred from other areas across inner time step
    # qin <- dist.flux(groups, 
    #                  flows$qbf,
    #                  ichan = ichan,
    #                  W=w)
    # 
    
   	flows$qin <- as.vector((flows$qbf*groups$area) %*% w)

  	# update storage deficits with the new input and output flows
  	# deficit reduced by upslope input and unsat drainage and
  	# increased by downslope flow in to other units
   	# existing excess storage is removed from deficit
  	SD.add <- (flows$qbf - flows$qin/groups$area - flows$uz)*dtt - stores$ex
		SD.add[ichan] <- 0
  	stores$sd <- stores$sd + SD.add

  	# -ve SD = excess surface storage
  	sat.ex <- pmin(stores$sd, 0)
 	
  	# Deficit excess, note -ve
  	stores$ex <- -sat.ex #stores$ex    # + flows$ex*dtt
  	
  	# 
  	stores$sd <- pmax(stores$sd, 0)

    # actual ae at this time step
    ae.dtt[inner,] <- flows$ae
    
    fc <- flows[ichan, ]
  	# total flux transferred into river plus rainfall minus evap
    Qriv <- Qriv + (fc$qin + (fc$rain-fc$ae) * achan)/ntt
    
    # increase in deficit
 #   sd.add <- (flows$qbf - flows$qin/groups$area -  flows$uz)*dtt  
    
  	# record intermediate unsaturated drainage
  	uz.inner[inner,] <- flows$uz  
  }
	
#	stores$sd <- stores$sd + SD.add
	
	# total drainage over time step (catches situation where unsat zone drains over entire period)
	flows$uz <- signif(colMeans(uz.inner), 3)
	
  # average out total ae
  flows$ae <- colMeans(ae.dtt)  #Sums(ae.dtt)*dtt/dt

  # total input into channel
  flows[ichan,]$qin <- Qriv
  
  # specific overland flow (rate)
#  flows$qof <- stores$ex/dt

	# ############################## water balance check ###########################
# 	stores.diff <- stores - stores0
# 	store.gain <- stores.diff$ex + stores.diff$srz + stores.diff$suz-stores.diff$sd
# 	net.in <- (flows$rain- flows$ae)*dtt
# 	bal <- store.gain-net.in
	# ############################## water balance check ###########################
  # return updated fluxes and storages, total discharge into river
  return(list("flows"=flows, "stores"=stores))
}

# adjust storages given updated base flow and inputs over this time step
# if max storage deficit or saturation recahed due to net inflow (outflow) then
# route the excess overland
update.storages <- function(groups, flows, stores, dtt, ichan, tm)
{
  # initially add any excess baseflow to the excess(surface) storage
  stores$ex <- stores$ex + flows$ex*dtt

  noflow <- setdiff(which(stores$sd>=groups$sd_max), ichan)
  if(length(noflow)>0)
  {
    LogEvent("SD > SDmax")  #, tm=tm)  #, paste0(groups[noflow,]$id, collapse=",")), tm=tm)
    stores[noflow,]$sd<- groups[noflow,]$sd_max
    #cat("Maximum storage deficit reached. This usually indicates pathological behaviour leading to extreme performance degradation. Execution halted")
   # stop()
    flows[noflow,]$qbf <- 0  # flows[noflow,]$qbf - (stores[noflow,]$sd - groups[noflow,]$sd_max) / dtt
  }
	# check for max saturated flow for grouping, in which case set SD to zero and
	# route excess as overland flow update storage deficits using the inflows and
	# outflows (fluxes) - note that sd is a deficit so is increased by base flow
	# out of groups and increased by flow from unsaturated zone and upslope areas
	# inc drainage from unsat zone
	# net outflow from land zone - adds to sd.
	bal <- dtt*(flows$qbf - flows$uz -flows$qin/groups$area)

	# inflow / outflow for channel is handled by the time delay histogram so storage isn't relevant here
	bal[ichan]<-0
	stores$sd <- stores$sd + bal


	# do not allow base flow to reduce storage  below zero. This should be routed overland
	saturated <- which(stores$sd < 0)  #    #setdiff(which(stores$sd<=0), ichan)
	if(length(saturated)>0)
	{
		#LogEvent(paste("Base flow saturation in zones ", paste(groups[saturated,]$tag, collapse=",")), tm=time)
		#  browser()
		# transfer negative deficit to excess store
		stores[saturated,]$ex <- stores[saturated,]$ex  -  stores[saturated,]$sd
		#
		#     # add to overland storage
		#     flows[saturated,]$ex <-
		stores[saturated,]$sd<- 0
	}

	# check channels
	return(list("flows"=flows, "stores"=stores))
}


