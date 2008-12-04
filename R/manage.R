
manage <- function(init, epistep, vacgrid, costs, T=40,
                   pinit = list(b=.1,k =.02,nu=.2,mu=.1),          # initial parameter values
                   hyper = list(bh=c(1,3), kh=c(1,3), nuh=c(1,1), muh=c(1,1)), #hyperparameter
                   vac0 = list(frac=0, stop=0),                    # policy enacted before estimation, if any
                   MCvits=10, MCMCpits=1000, vacsamps=100,
                   start=8 , ...)                              # time step to begin vaccination calculation
{
  ##   graphics.off();
  ##   x11(); x11();
  ##   mydevs = dev.list()
  howmany=0
  doagain=TRUE
  while(doagain && howmany<100)
    {      
      ## initialize the soln data frame
      soln <- init.soln(init, T)
      
      ## initialize the prime data frame
      prime <- data.frame(matrix(0, ncol=2, nrow=T))
      names(prime)<-c("S", "I")
      
      ## initialize the vaccination and culling policies
      VAC<- vac0$frac
      STOP<-vac0$stop
      allVfracs <- NULL  ## going to store the distribution of vfracs across samples of parameter values
      allVstops <- NULL  ## going to store the distribution of vstops across samples of parameter values
      allVtimes <- NULL
      allpolicies <- matrix(c(0,0), nrow=1)
      havesamps=FALSE
      samp<-NULL
      ## starting values for the SIR parameters
      b <- pinit$b; k <- pinit$k; nu <- pinit$nu; mu <- pinit$mu;
      ##   b<-1; k<-1
      ##   nu <- 0.1; mu <- 0.001
      
      doagain=FALSE
      ## step through time
      for( i in 2:T ){
        
        ## if you get to the start time, and the epidemic hasn't grown, start over
        if(i==start && soln$I[i-1]<=soln$I[1]) 
          {
            doagain=TRUE
            break;
            
          }
        
        ## check if there is a need to vaccinate
        if(! (soln$S[i-1] == 0 || soln$I[i-1] == 0) && i>=start && havesamps && soln$S[i-1]>STOP)
          {
            ## call newvacpolicy so that vcsamps samples are taken
        ## conditional on a thinned version of samp
            VACs <- vaccinate(soln$S[i-1], soln$I[i-1], samp, vacsamps,
                               costs, MCvits, T-i+1, vacgrid)
            
            ##       ## plot a histogram of the vaccination policy
            ##       dev.set(mydevs[2])
            ##       ## add a plot for the cost surface here.
            ##       image(Vprobs,Vstops,VACs$bestC,xlab="fraction",
            ##             ylab="stop number",main="Best Policy Histogram");
            ##       grid(length(Vprobs),length(Vstops),lty=1,col='black');box()
            
            ## concatenate the policy for this time step to the history
            allVfracs <- rbind(allVfracs, VACs$Vfractions)
            allVstops <- rbind(allVstops, VACs$Vstoptimes)
            allVtimes <- c(allVtimes,i)
            ##      Vold <- VACs$Vold
            
            ## make the current vaccination policy equal to the mean
            mostbestpol = which(VACs$bestC==max(VACs$bestC),arr.ind=TRUE)[1,]
            VAC <- vacgrid$fracs[mostbestpol[1]]##mean(VACs$Vfractions)
            STOP<- vacgrid$stops[mostbestpol[2]]##mean(VACs$Vstoptimes)
          }
        
        ## run epidemic one step forward
        allpolicies<-rbind(allpolicies, c(VAC,STOP))
        
        ## simulate one step forward in the epidemic, with epistep, with
        ## vaccinations and cullings coming first in epimanage
        out <- epimanage(soln=soln, epistep=epistep, i=i, VAC=VAC, STOP=STOP, ...)
        
        ## update (prime) totals
        prime[i,] <- epi.update.prime(soln[i-1,], out)
        
        ## update (soln) totals
        soln[i,] <- epi.update.soln(soln[i-1,], out, costs)
        
        ## do T MH/Gibbs steps to sample from the join distrib of b, k, mu and nu
        if(soln$S[i]>0 && soln$I[i]>0)
          {
            samp <- mcmc.bknumu(MCMCpits, b, k, nu, mu, soln$itilde[2:i], soln$rtilde[2:i],
                                soln$dtilde[2:i], prime$S[2:i], prime$I[2:i], hyper)
            
            ## collect means    
            b <- mean(samp$b); k <- mean(samp$k); nu <- mean(samp$nu);
            mu <- mean(samp$mu)
            havesamps=TRUE
          }
      }
    }	
  ##X11()
  ##png(filename="novac_solnplot.png")
  ##plotsolns(soln)
  ##dev.off()
  
  ##  x11();TimeSeriesOfDensities(allVfracs,"vaccination fraction")
  ##  x11();TimeSeriesOfDensities(allVstops,"stop number")
  ## soln and the (full) calculated policies
  ##  x11();PlotCosts(soln,main="Adaptive Vaccination");box()

  ## make allpolicies into a data frame
  allpolicies <- data.frame(allpolicies)
  names(allpolicies) <- c("frac", "stop")
  
  r <- list(soln=soln, vachist=list(fracs=allVfracs, stops=allVstops),
            pols=allpolicies, vactimes=allVtimes, samp=samp)
  r$call <- match.call()
  class(r) <- "epiman"
  return(r)
}


## init.soln:
##
## function to initialize the solution -- the data frame
## that contains all of the information about the epidemic
## -- starting with a particular number of succeptables and
## infecteds, for a total time horizon t

init.soln <- function(init, T)
  {
    soln <- data.frame(matrix(0, ncol=11, nrow=T))
    names(soln)<-c("TIME", "S", "I", "R", "D", "itilde", "rtilde", "dtilde",
                   "V", "QC","C")
    soln$TIME[1] <- 1        # column 1: time
    soln$S[1] <- init$S     # column 2: susceptibles
    soln$I[1] <- init$I     # column 3: infecteds

    return(soln)
  }


## epi.update.prime:
##
## the real number of succeptiables and infecteds that are
## around at the beginning of the current time step, where
## soln contains the state at the previous timestep and
## out contains information about the number just culled
## and just vaccinated

epi.update.prime <- function(soln, out)
  {
    prime <- data.frame(S=soln$S - out$justvacc,
                        I=soln$I - out$justcull)
    return(prime)
  }


## epi.update.soln:
##
## advance time by one unit in the soln data frame based
## on the previous time's state information (contained in s)
## and information (in out) about the number just culled,
## vaccinated, died, recovered, newly infected, etc

epi.update.soln <- function(s, out, costs)
  {
    soln <-
      data.frame(
                 TIME = s$TIME+1,
                 S = s$S - out$justvacc - out$newi,
                 I = s$I - out$justcull - out$justrec - out$justdied+out$newi,
                 R = s$R + out$justrec, #+ out$justvacc + out$justcull
                 D = s$D + out$justdied,
                 itilde = out$newi,
                 rtilde = out$justrec,
                 dtilde = out$justdied,
                 V = s$V + out$justvacc,  
                 QC = s$QC + out$justcull,
                 C = s$C+ costs$infect*s$I+costs$vac*out$justvacc+costs$death*out$justdied)
    
    return(soln)
  }


## vaccinate:
##
## function to call newvacpolicy on a thinned version of the
## samples contained in samp.  The thinning is determined by
## the desired number of samples from the vaccination policy
## (vacsamps).  nI and nS are the number of infecteds and
## succeptibles in the current timestep.  Vold is the last
## used vaccination policy.  the rest of the arguments are
## identical to the manage function.

## so that I can step through the vaccinate function
## nS=soln$S[i];nI=soln$I[i];samp=samp;vacsamps=vacsamps;
## ci=ci;cv=cv;cd=cd;MCvits=MCvits;timehorizon=t-i+1;Vprobs=Vprobs;Vstops=Vstops;

vaccinate <- function(nS,nI, samp, vacsamps, costs,
                      MCvits, timehorizon, vacgrid)
  {
    ## cat("finding vaccination strategy for (S,I)=(",nS,",",nI,"),
    ##     ",timehorizon," days left\n",sep="")
    
    ## total number of active individuals
    N <- nI + nS

    ## indices of the thinned sample
    ll<-floor(seq(1,nrow(samp), length.out=vacsamps))#by=floor(nrow(samp)/vacsamps))

    ## the thinned sample
    use.samp<-samp[ll,]
    Vfractions<-rep(0,length(ll))
    Vstoptimes<-rep(0,length(ll))
    Vbestcosts = matrix(0,length(vacgrid$fracs),length(vacgrid$stops))
    Vcosts = array(0,c(length(ll),length(vacgrid$fracs),length(vacgrid$stops)))
    ## get the vaccination policy for each thinned sample,
    ## starting the next call where we left off on the
    ## last call
    for(j in 1:length(ll)){

      ## call the function which gets the sample
      VP <- VarStopTimePolicy(nS,nI,timehorizon,use.samp$b[j],
                              use.samp$k[j],use.samp$nu[j],use.samp$mu[j],
                              costs$vac, costs$death, costs$infect,
                              MCvits, vacgrid$frac, vacgrid$stops,
                              midepidemic=TRUE,start=0)
      bestpol = which(VP==min(VP),arr.ind=TRUE)[1,]
      Vfractions[j]=vacgrid$fracs[bestpol[1]]
      Vstoptimes[j]=vacgrid$stops[bestpol[2]]
      Vbestcosts[bestpol[1],bestpol[2]] = Vbestcosts[bestpol[1],bestpol[2]]+1
      Vcosts[j,,] = VP
      ## cat("  j=",j," b=",use.samp$b[j]," k=",use.samp$k[j]," nu=",use.samp$nu[j],
      ##    " mu=",use.samp$mu[j]," vf=",Vfractions[j]," vs=",Vstoptimes[j],"\n",sep="")
      if(max(Vbestcosts)>vacsamps/2)
        {
          break
        }
      ##       VP<-newvacpolicy(N=N, b=use.samp$b[j], k=use.samp$k[j],
      ##                        mu=use.samp$mu[j], nu=use.samp$nu[j],
      ##                        q=q, ci=ci, cv=cv, cd=cd, MCvits=MCvits,
      ##                        V.in=Vold, thresh=thresh)
      ##       MCvits[1]<-max(100,VP$MCvits[1]/10)
      ##       Vold<-VP$V
      
      ## get the sample
      
      ##cat(paste("j=", j, " iters=", VP$iters, "\n", sep=""))
      ##       VACs[j]<-VP$alpha[nI+1,nS+1]
    }

    return(list(Vfractions=Vfractions,Vstoptimes=Vstoptimes,C=Vcosts,bestC=Vbestcosts))
  }


## epistep:
##
## simulate the epidemic one time step foreward based on
## the previous time step and the parameters k, b, nu and mu

epistep <- function(SIR, true=list(b=0.00218, k=10, nu=0.4, mu=0))
  {
    b <- true$b; k <- true$k; nu <- true$nu; mu <- true$mu
    
    ## compute the number of removed
    rem <- rbinom(1, SIR$I, 1-exp(-(nu + mu)))
    
    ## compute the number of recovered
    rec <- rbinom(1, rem, nu/(nu + mu))
    
    ## compute the binomial probability of a new infection
    ## can use soln[i-1,3] or soln[i-1,3]-justcull
    p <- 1-(k/(k+b*(SIR$I)))^k   
    
    ## sample from the binomial distribution
    infect <- rbinom(1, SIR$S, p)
    
    ## compute the number which just died
    dead <- rem - rec
    
    return(list(rem=rem, rec=rec, infect=infect, dead=dead))
  }


## epimanage:
##
## simulate the epidemic one time step forward.  After managing
## via vaccinations and culling, then apply epistep on
## whatever is left.  The elipses argument (...) should contain
## whatever epistep needs to do its business, like b,k,nu, and mu
## (true) for example

epimanage <- function(soln, epistep=epistep, i, VAC=0.1, STOP=1, ...)
{ 
    ## compute the vacc/cull policy
    ##  if( i == 2){  ## why is this set at 2? maybe we want to vaccinate immediately
    ##     	vacc<-0
    ##     	cull<-0	
    ##     }else{
    if(soln$S[i-1]>STOP) vac <- round(VAC*soln$S[i-1])
    else vac <- 0
    ##  	justcull<-round(CULL*soln$I[i-1])
    cull<-0
    ##}
    
    ## use epistep to simulate one time step forward based on the
    ## previous time step and the parameters b, k, nu, and mu
    SIR <- list(S=soln$S[i-1]-vac, I=soln$I[i-1]-cull, R=NA)
    just <- epistep(SIR, ...)

    ## assemble the outputs
    out <- matrix(0, ncol=5, nrow=1)
    out[1,1] <- vac
    out[1,2] <- cull
    out[1,3] <- just$rec
    out[1,4] <- just$infect
    out[1,5] <- just$dead
    out <- data.frame(out)
    names(out) <- c("justvacc", "justcull", "justrec", "newi", "justdied")
    out
  
}
