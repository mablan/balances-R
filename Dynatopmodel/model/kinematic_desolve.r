# downslope distribution function resulting from exponential transmissivity
# assumption
dqdt <- function(q, times, w, m, r, a)
{
	q <- y
#	dddt <- q- (t(w)%*%(a*q))/a-r
#	dqdt <- dddt*-q/m
	A <- diag(1/a) %*% t(w) %*% diag(a) -identity.matrix(nrow(w))
  # throttle input
  #q <- pmax(q, qmax)
	res <- q/m * (A %*%q + r)
	return(res)
#	return(res)
}

dqdt.ode <- function(t, y, parms, 
                      ...)
{
  if(length(parms$idry)>0)
  {
  	y[parms$idry] <- 0	
  }
  
	#	dddt <- q- (t(w)%*%(a*q))/a-r
	#	dqdt <- dddt*-q/m
#	A <- diag(1/a) %*% t(w) %*% diag(a) -identity.matrix(nrow(w))
	# 

	res <-  y/parms$m * (parms$A %*% y + parms$r)
		
# 	iex <- which((res*dt+y) > parms$qbmax)
# 	#	-fun(q, S, params=parms$groups) * (A %*%q + r)
#   # impose maximum q and ensure non-negative
# 	if(length(iex)>0)
# 	{
# 	  res[iex] <- 0
# 	}

# 	if(any(res > parms$qbmax, na.rm=TRUE))
# 	{
#   res <- pmin(res, parms$qbmax, na.rm=TRUE)  #pmax(pmin(res, parms$groups$qbmax, na.rm=T),0)
# 	}
  # ex <- pmax(res-res2, 0)
	return(list(res))
	#	return(res)
}

#cdqdt.ode <- cmpfun(dqdt.ode)

require(compiler)
require(deSolve)

# ******************************************************************************
# kinematic.r kinematic wave routes flow downslope through groupings from
# 4-point solution using input flows and output flows from previous and current
# time steps
# ==============================================================================
# Inputs
# ------------------------------------------------------------------------------
# time step dt in hrs
# flowst1: group flows and storages at previous time step
# flows: groups flows at current time step- Qbf to be determined
# w: time-stepping weighting coefficient: zero for a totally implicit solution
# (depends only on flows for previous steps). w=0.5
# niter: max number of iterations in iterative scheme
# ==============================================================================
# Returns
# ------------------------------------------------------------------------------
# updated storage deficts, estimates for base flows at this time step
# ==============================================================================
# References
# ------------------------------------------------------------------------------
# Beven and Freer (2001). A Dynamic TOPMODEL
# Beven (1981). Kinematic subsurface stormflow
# Li el al (1975).  Li, Simons and Stevens 1975 - Nonlinear Kinematic Wave
# Approximation for Water Routing
# Beven (2012). Rainfall runoff modelling, chapter 5 pp141-150, pp.180-183
# ******************************************************************************
# note: we now exclude lateral input from land hsus from the flux inpu
route.kinematic.euler <- function(groups,
                            flows,        # fluxes at previous time step (prediction time)
                            stores,       # current storage
                            dtt,
                            ichan,
                            time,
							              w,
                            nstep=1,
				                    dqds,
														A=NULL,
                            method="lsoda"              # Livermore solver
)
{
  # UZ: recharge is drainage from unsaturated into saturated zone - assumed constant over time steps
  r <- flows$uz #

 	qb0 <- flows$qbf
 	qb0[ichan] <- 0
 	
 	
 	#	dddt <- q- (t(w)%*%(a*q))/a-r
 	#	dqdt <- dddt*-q/m
 	if(is.null(A))
 	{
 		a <-groups$area
 		#  fun <- parms$dqds
 		
 		N <- nrow(w)
 		
 		 A <- diag(1/a, N, N) %*% t(w) %*% diag(a, N, N) - identity.matrix(N)
 	}
 	
 	# identify any units that have dried out andare not producing any base flow
 	idry <- which(stores$sd>=groups$sd_max)
 	
 	# nstep is the number of results to produce
 	res <- ode(y=qb0, 
 						 times=seq(0, dtt, length.out=nstep+1),
             func=dqdt.ode,
 						 parms=list(A=A,
 						 						idry=idry,
 						 					  m=groups$m,
 						 						qbmax=groups$qbmax,
 						 						r=r),
 						 method=method)
 	
  # final row gives baseflows at the end stage
 	flows$qbf <- res[nstep+1,-1]

#  res <- matrix(res[-1,-1], nrow=nstep)
  
 # flows$ex <- 0# pmax(flows$qbf-groups$qbf_max, 0)
  
 	iex <- which(flows$qbf > groups$qbmax)
   if(length(iex) >0)
   {
     flows$qbf[iex] <- groups$qbmax[iex]
   }
	# throttle outlet discharge to max given gradient and conductivity and ensure >0
 #	flows$qbf <- flows$qbf-flows$ex
 	
 	# flow out of areas that have reached their maximum deficit are zero
 	
 	
  return(flows)
}





