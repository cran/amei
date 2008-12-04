SimulateEpidemic = function(S0,I0,R0,D0,T,b,k,nu,mu,vacc,vaccstop,cvacc,
                            cdeath,cinfected,starttime)
  {
    cout = .C("RSimulateEpidemic",
      S = as.integer(rep(S0,T)),
      I = as.integer(rep(I0,T)),
      R = as.integer(rep(R0,T)),
      D = as.integer(rep(D0,T)),
      V = as.integer(rep(0,T)),
      C = as.double(rep(0,T)),
      T = as.integer(T),
      b = as.double(b),
      k = as.double(k),
      nu = as.double(nu),
      mu = as.double(mu),
      vacc=as.double(vacc),
      vaccstop=as.integer(vaccstop),
      cvacc = as.double(cvacc),
      cdeath = as.double(cdeath),
      cinfected = as.double(cinfected),
      starttime=as.integer(starttime),
      PACKAGE="amei")

    list(S = cout$S,I=cout$I,R=cout$R,D=cout$D,V=cout$V,C=cout$C)
  }

## ok, if midepidemic is true, then we're just going to accept all epidemics
## if midepidemic is false, then we're going to assume we're at the beginning 
## of an epidemic, and that starttime is something >0.  in this case, we'll 
## throw away any epidemics for which the number of infecteds hasn't grown 
## by startime. 
VarStopTimePolicy = function(S0,I0,T,b,k,nu,mu,cvacc,cdeath,cinfected,
                             MCvits,Vprobs,Vstops,midepidemic,starttime)
  {
    cout = .C("RVarStopTimePolicy",
      S0 = as.integer(S0),
      I0 = as.integer(I0),
      T = as.integer(T),
      b = as.double(b),
      k = as.double(k),
      nu = as.double(nu),
      mu = as.double(mu),
      cvacc = as.double(cvacc),
      cdeath = as.double(cdeath),
      cinfected = as.double(cinfected),
      MCvits = as.integer(MCvits),
      Vprobs = as.double(Vprobs),
      nVprobs = as.integer(length(Vprobs)),
      Vstops = as.integer(Vstops),
      nVstops = as.integer(length(Vstops)),
      EC = double(length(Vprobs)*length(Vstops)),
      midepidemic=as.integer(midepidemic),
      starttime=as.integer(starttime),
      PACKAGE="amei")
    C <- matrix(cout$EC,length(Vprobs),length(Vstops),byrow=TRUE)
    return(C)
  }


SimulateManagementQuantiles <- function(epistep,Time,init, pinit, hyper, vac0,
                                        costs, start, MCvits, MCMCpits,vacsamps,
                                        vacgrid, nreps,lowerq, upperq, verb=FALSE, ...)
  {
    Sall = matrix(0,nrow=Time,ncol=nreps)
    Iall = matrix(0,nrow=Time,ncol=nreps)
    Rall = matrix(0,nrow=Time,ncol=nreps)
    Dall = matrix(0,nrow=Time,ncol=nreps)
    Vall = matrix(0,nrow=Time,ncol=nreps)
    Call = matrix(0,Time,nreps)
    PoliciesAll = array(0,c(Time,2,nreps))
    
    for(n in 1:nreps)
      {
        if(verb) cat("*** Simulating epidemic",n,"***\n")
        foo <- manage(epistep=epistep,pinit=pinit,T=Time,init=init,hyper=hyper,
                          vac0, costs=costs, MCMCpits=MCMCpits,
                          vacsamps=vacsamps, vacgrid=vacgrid, start=start, ...)
        Sall[,n] = foo$soln$S
        Iall[,n] = foo$soln$I
        Rall[,n] = foo$soln$R
        Dall[,n] = foo$soln$D
        Vall[,n] = foo$soln$V
        Call[,n] = foo$soln$C
        PoliciesAll[,,n] = as.matrix(foo$pols)
      }
    
    SQ1 = apply(Sall,1,quantile,prob=lowerq)
    Smean = apply(Sall,1,mean)
    Smed = apply(Sall,1,median)
    SQ3 = apply(Sall,1,quantile,prob=upperq)
    
    IQ1 = apply(Iall,1,quantile,prob=lowerq)
    Imean = apply(Iall,1,mean)
    Imed = apply(Iall,1,median)
    IQ3 = apply(Iall,1,quantile,prob=upperq)
    
    RQ1 = apply(Rall,1,quantile,prob=lowerq)
    Rmean = apply(Rall,1,mean)
    Rmed = apply(Rall,1,median)
    RQ3 = apply(Rall,1,quantile,prob=upperq)
    
    DQ1 = apply(Dall,1,quantile,prob=lowerq)
    Dmean = apply(Dall,1,mean)
    Dmed = apply(Dall,1,median)
    DQ3 = apply(Dall,1,quantile,prob=upperq)
    
    VQ1 = apply(Vall,1,quantile,prob=lowerq)
    Vmean = apply(Vall,1,mean)
    Vmed = apply(Vall,1,median)
    VQ3 = apply(Vall,1,quantile,prob=upperq)
    
    CQ1 = apply(Call,1,quantile,prob=lowerq)
    Cmean = apply(Call,1,mean)
    Cmed = apply(Call,1,median)
    CQ3 = apply(Call,1,quantile,prob=upperq)
    
    PolQ1 = apply(PoliciesAll,c(1,2),quantile,prob=lowerq)
    Polmean = apply(PoliciesAll,c(1,2),mean)
    Polmed = apply(PoliciesAll,c(1,2),median)
    PolQ3 = apply(PoliciesAll,c(1,2),quantile,prob=upperq)
    
    list(Q1 = data.frame(S=SQ1,I=IQ1,R=RQ1,D=DQ1,V=VQ1,C=CQ1,frac=PolQ1[,1],stop=PolQ1[,2]),
         Mean = data.frame(S=Smean,I=Imean,R=Rmean,D=Dmean,V=Vmean,C=Cmean,frac=Polmean[,1],stop=Polmean[,2]),
         Median = data.frame(S=Smed,I=Imed,R=Rmed,D=Dmed,V=Vmed,C=Cmed,frac=Polmed[,1],stop=Polmed[,2]),
         Q3 = data.frame(S=SQ3,I=IQ3,R=RQ3,D=DQ3,V=VQ3,C=CQ3,frac=PolQ3[,1],stop=PolQ3[,2]))       
  }


  

SimulateEpidemicQuantiles = function(S0,I0,R0,D0,T,b,k,nu,mu,vacc,vaccstop,cvacc,cdeath,
  cinfected,nreps,lowerq,upperq,midepidemic,starttime)
  {
    Sall = matrix(0,nrow=T,ncol=nreps)
    Iall = matrix(0,nrow=T,ncol=nreps)
    Rall = matrix(0,nrow=T,ncol=nreps)
    Dall = matrix(0,nrow=T,ncol=nreps)
    Vall = matrix(0,nrow=T,ncol=nreps)
    Call = matrix(0,T,nreps)
    for(n in 1:nreps)
      {
        blah=TRUE;blahcount=0
        while(blah)
          {
            foo = SimulateEpidemic(S0,I0,R0,D0,T,b,k,nu,mu,vacc,vaccstop,cvacc,cdeath,cinfected,starttime)
            if(!midepidemic)
              {
                blahcount = blahcount+1
                if(foo$I[starttime-1]>I0)##sum(foo$I[-1]>foo$I[1])>1 | (b<=0 | k<=0) | vacc==1 | midepidemic)
                  blah=FALSE
                if(blahcount==100)
                  {
                    cat("Warning: <1% chance of an epidemic\n")
                    blah=FALSE
                  }
              }
                
          }
        Sall[,n] = foo$S
        Iall[,n] = foo$I
        Rall[,n] = foo$R
        Dall[,n] = foo$D
        Vall[,n] = foo$V
        Call[,n] = foo$C
      }

    SQ1 = apply(Sall,1,quantile,prob=lowerq)
    Smean = apply(Sall,1,mean)
    Smed = apply(Sall,1,median)
    SQ3 = apply(Sall,1,quantile,prob=upperq)

    IQ1 = apply(Iall,1,quantile,prob=lowerq)
    Imean = apply(Iall,1,mean)
    Imed = apply(Iall,1,median)
    IQ3 = apply(Iall,1,quantile,prob=upperq)

    RQ1 = apply(Rall,1,quantile,prob=lowerq)
    Rmean = apply(Rall,1,mean)
    Rmed = apply(Rall,1,median)
    RQ3 = apply(Rall,1,quantile,prob=upperq)

    DQ1 = apply(Dall,1,quantile,prob=lowerq)
    Dmean = apply(Dall,1,mean)
    Dmed = apply(Dall,1,median)
    DQ3 = apply(Dall,1,quantile,prob=upperq)
    
    VQ1 = apply(Vall,1,quantile,prob=lowerq)
    Vmean = apply(Vall,1,mean)
    Vmed = apply(Vall,1,median)
    VQ3 = apply(Vall,1,quantile,prob=upperq)

    CQ1 = apply(Call,1,quantile,prob=lowerq)
    Cmean = apply(Call,1,mean)
    Cmed = apply(Call,1,median)
    CQ3 = apply(Call,1,quantile,prob=upperq)

    
    list(Q1 = data.frame(S=SQ1,I=IQ1,R=RQ1,D=DQ1,V=VQ1,C=CQ1),
         Mean = data.frame(S=Smean,I=Imean,R=Rmean,D=Dmean,V=Vmean,C=Cmean),
         Median = data.frame(S=Smed,I=Imed,R=Rmed,D=Dmed,V=Vmed,C=Cmed),
         Q3 = data.frame(S=SQ3,I=IQ3,R=RQ3,D=DQ3,V=VQ3,C=CQ3))       
  }
