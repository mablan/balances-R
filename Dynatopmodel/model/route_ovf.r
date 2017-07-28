dsdt <- function(t, st, parms, ...)
{
  # loss of specific storage per uni time (each unit)
  qout <-  parms$groups$vof*st
  
  # gain of specific storage per unit time
  
  qin <- ((qout*parms$groups$area) %*% parms$W)/parms$groups$area
  
  dsdt <- qin - qout

  return(list(dsdt))
  
}


# routing of surface excess  using a system of ODEs constructed with flux routing matrix
distribute_surface_excess_storage <- function(groups, 
                           W, 
                           ex, 
                           dt, 
                           ichan=1)
{
  if(any(ex > 0, na.rm=T))
  {

    # vof a vector of overflow flow
    #vof <- groups$vof
    #groups$vof[ichan] <- 0
    
    # relate storage to discharge by q = vs
    # then s' = t(q) %*% w-q = v*(t(w)%*% s-s)
    #	ex.hru[ichan]<- 0
    #		eig <- get_routing_eig(A, vof)
    
    # A is the scaled weighting matrix for this discretisation
    ex.dtt <- ode(y=ex,
                  times=c(0, dt), func=dsdt, 
                  parms=list(W=W, groups=groups))
    
    ex <- round(ex.dtt[2,-1],5)
    
    # anything that can't be accommodated on surface
    ex.over <- pmax(ex[]-groups$ex_max, 0, na.rm=TRUE)
    
    # send overflow directly to channel
    if(any(ex.over[-ichan]>0))
    {
      ex[ichan] <- ex[ichan]+ sum(ex.over*groups$area)/groups$area[ichan]
      ex <- ex - ex.over
    }
  }
  return(ex)
}

simple_ovf_routing <- function(groups, ex, W, dt, ichan=1)
{
#  ex.hru <- round(ex, )  # ex[-ichan]
  if(any(ex > 0, na.rm=T))
  {
    # total 
    Ex <- ex*groups$area

    # total flux out of each group
    Qout <- groups$vof*Ex
    
    # distributed to downslope groups
    Qin <- Qout %*% W 
    
    # change in storage over time step
    Ex <- pmax(Ex + (Qin - Qout)*dt, 0)
    
    ex <- as.vector(Ex/groups$area)
  }
  return(ex)
    
}

# use the eigenvector approach to solving a system of linear first order ODEs
dist.eigen <- function(groups, A, ex, dt, ichan=1, eig=NULL, tm=NULL, max.t=NULL)
{

	if(any(ex > 1e-6, na.rm=T))
	{
		ex.hru <- ex #* groups$area
		
		# vof a vector of overflow flow
		vof <- groups$vof
		
		# ensure no storage routed out of channel
		vof[ichan] <- 0
		
		# relate storage to discharge by q = vs
		# then s' = t(q) %*% w-q = v*(t(w)%*% s-s)
	#	ex.hru[ichan]<- 0
#		eig <- get_routing_eig(A, vof)
		
		# A is the scaled weighting matrix for this discretisation
		ex.dtt <- ode(y=ex.hru,
                  times=c(0, dt), func=dsdt, 
		              parms=list(A=A, vof=vof))
		
		ex.hru <- ex.dtt[2,-1]  # exclude first column with times and take result at time dt, in final row

#		Av <- vof * (A - identity.matrix(length(vof)))
#		eig <- eigen(Av)

		# scale to give rate of change of storage in tewrms of storgage

	#	vof[ichan] <- 0 #vof[-ichan]
#		vof <- matrix(rep(vof), col=nrow(w))
# 		A <- t(whru) -identity.matrix(nrow(whru))
#
# 		A <- vof*A
# 		eig <- eigen(A)
		# eigenvalues
#		lambda <- eig$values
#		u <- eig$vectors

		# the generral asolution to this system is
		# s = c1*exp(l1*t)u + c2*exp(l2*t*u)+...

		# where the ci are arbitrary constants


		# Assume that ex is the initial state of the system so, setting t=0
		# ex = c[i]u[i]=t(c)%*%u=t(u)%c

		# solve for coefficients using initial conditions
#		ci <- solve(u, ex)

		# factors at t=dt
#		cdt <- ci * exp(lambda*dt)

#		ex.hru <- u %*% cdt

		# convert back to specific storage

		# the channel storages are unchanged
#		ex[-ichan] <- ex.hru
# specific storages
  #  ex.hru<- 0
	  ex <- ex.hru
	}

	# convert storage back to equivalent flow on return
	return(ex)
}

# dsdt <- function(t, st, parms, ...)
# {
# 	A <- parms$A
# 	vof <- parms$vof
#     # convert to a discharge
# 	Av <- vof * (A - identity.matrix(length(vof)))
#   dsdt <- vof*st %*% t(A)-vof*st
#   # Av %*% s
# 	return(list(dsdt))
# 
# }




# method of solution by eigen values and vectors
get_routing_eig <- function(A, vof, ichan=1)
{

	# scale to give rate of change of storage in terms of storage
	nhru <- nrow(A)

	#	vof[ichan] <- 0 #vof[-ichan]
	vof.mat <- matrix(rep(vof, nhru), nrow=nhru, byrow=T)
	Av <- A - diag(vof)

	eig <- eigen(Av)

	return(eig)

}
