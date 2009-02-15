###################################################
### chunk number 1: 
###################################################
library(amei)
options(width=60)


###################################################
### chunk number 2: 
###################################################
seed <- 12345
set.seed(seed)


###################################################
### chunk number 3: 
###################################################
true <- list(b=0.00218, k=10, nu=0.4, mu=0) 
init <- list(S0=762, I0=1, R0=0, D0=0) 
costs <- list(vac=2, death=4, infect=1)


###################################################
### chunk number 4: 
###################################################
vac <- list(frac=0, stop=0)


###################################################
### chunk number 5: 
###################################################
init.MCepi <- MCepi(init, true, vac, costs)


###################################################
### chunk number 6: epis
###################################################
plot(init.MCepi)


###################################################
### chunk number 7: costs
###################################################
plot(init.MCepi, type="costs")


###################################################
### chunk number 8: 
###################################################
vacgrid <- list(fracs=seq(0,1.0,0.1), stops=seq(2,init$S0-75,75)) 


###################################################
### chunk number 9: 
###################################################
out.optvac <- optvac(init, true, vacgrid, costs) 


###################################################
### chunk number 10: 
###################################################
best <- getpolicy(out.optvac) 
worst <-getpolicy(out.optvac, which="worst") 
rbind(best, worst) 


###################################################
### chunk number 11: optvac
###################################################
plot(out.optvac) 


###################################################
### chunk number 12: 
###################################################
vac.opt <- best[3:4]
opt.MCepi <- MCepi(init, true, vac.opt, costs)


###################################################
### chunk number 13: episov
###################################################
plot(opt.MCepi)


###################################################
### chunk number 14: costsov
###################################################
plot(opt.MCepi, type="costs")


###################################################
### chunk number 15: 
###################################################
getvac(opt.MCepi)


###################################################
### chunk number 16: 
###################################################
T <- length(opt.MCepi$Median$C)
optC <- getcost(opt.MCepi)
initC <- getcost(init.MCepi)
data.frame(rbind(initC,optC), row.names=c("init", "opt"))


###################################################
### chunk number 17: 
###################################################
opt.MCepi


###################################################
### chunk number 18: 
###################################################
out.man <- manage(init, epistep, vacgrid, costs) 


###################################################
### chunk number 19: epi
###################################################
plot(out.man)


###################################################
### chunk number 20: cost
###################################################
plot(out.man, type="cost")


###################################################
### chunk number 21: 
###################################################
getcost(out.man)


###################################################
### chunk number 22: params
###################################################
true <- as.list(formals(epistep)$true)
plot(out.man, type="params",true=true) 


###################################################
### chunk number 23: 
###################################################
out.man


###################################################
### chunk number 24: 
###################################################
out.MCmanage <- MCmanage(init, epistep, vacgrid, costs)


###################################################
### chunk number 25: MCmanepis
###################################################
plot(out.MCmanage)


###################################################
### chunk number 26: MCmancosts
###################################################
plot(out.MCmanage, type="costs")


###################################################
### chunk number 27: 
###################################################
getvac(out.MCmanage)


###################################################
### chunk number 28: MCmanfracs
###################################################
plot(out.MCmanage, type="fracs")


###################################################
### chunk number 29: MCmanstops
###################################################
plot(out.MCmanage, type="stops")


###################################################
### chunk number 30: 
###################################################
cinit <- getcost(init.MCepi)
copt <- getcost(opt.MCepi)
cman <- getcost(out.MCmanage)
data.frame(rbind(cinit, copt, cman), 
           row.names=c("init", "opt", "man"))


###################################################
### chunk number 31: 
###################################################
bad <- list(b=0.001, k=10, nu=0.9, mu=0)


###################################################
### chunk number 32: 
###################################################
costs.bad <- optvac(init, bad, vacgrid, costs)
pol.bad <- getpolicy(costs.bad)
pol.bad


###################################################
### chunk number 33: 
###################################################
bad.MCepi <- MCepi(init, true, pol.bad[3:4], costs)
cbad <- getcost(bad.MCepi)
cbad


###################################################
### chunk number 34: 
###################################################
seed <- 12345
set.seed(seed)


###################################################
### chunk number 35: 
###################################################
alt.epistep <- 
function(SIR, last=list(rem=0, rec=0, infect=0, dead=0, Z=0),
         true=list(a = 0.05, mu = 0.05, nu = 0.1, m = 0.4, 
	 rho = 200, C = 500))
{
  ## calculate the infection probability based on the
  ## resevoir, and randomly infect susceptibles
  Z <- last$Z
  fz <- Z/(Z+true$C)
  pi <- 1 - exp(-true$a * fz)
  infect <- rbinom(1, SIR$S, pi)

  ## update recovereds and deaths
  pr <- 1 - exp(-true$nu)
  rec <- rbinom(1,SIR$I,pr)
  pd <- 1 - exp(-true$mu)
  dead <- rbinom(1, SIR$I-rec, pd)

  ## resevoir dynamics
  pz <- 1 - exp(-true$m)
  dz <- rbinom(1, Z, pz)
  bz <- round(SIR$I*true$rho)
  Z <- Z - dz + bz

  ## the returned list is passed in as "last" in a
  ## subsequent call to this "epistep" function
  return(list(rem=(rec+dead), rec=rec, infect=infect, 
              dead=dead, Z=Z))
}


###################################################
### chunk number 36: 
###################################################
init1 <- list(S0=150, I0=1, R0=0, D0=0)
true<- list(a=0.065, mu=0.0, nu=0.3, m=0.99, rho=500, C=500)
alt.epistep1 <- alt.epistep
formals(alt.epistep1)$true <- true
out.alt<- manage(init1, alt.epistep1, NULL, NULL, T=80)


###################################################
### chunk number 37: 
###################################################
true <- list(a=0.065, mu=0.0, nu=0.3, m=0.005, rho=500, C=500)
alt.epistep2 <- alt.epistep
formals(alt.epistep2)$true <- true
out.alt2 <- manage(init1, alt.epistep2, NULL, NULL, T=80)


###################################################
### chunk number 38: alt
###################################################
plot(out.alt)


###################################################
### chunk number 39: alt2
###################################################
plot(out.alt2)


###################################################
### chunk number 40: alt-params
###################################################
plot(out.alt, type="params", showd=TRUE)


###################################################
### chunk number 41: alt2-params
###################################################
plot(out.alt2, type="params", showd=TRUE)


###################################################
### chunk number 42: 
###################################################
init <- list(S0=600, I0=1, R0=0, D0=0)
time=80
posterior <- manage(init, alt.epistep, NULL, NULL, T=time, bkrate=100)


###################################################
### chunk number 43: alt-mixing-b
###################################################
plot(log(posterior$samp$b), type="l", main="")


###################################################
### chunk number 44: alt-mixing-k
###################################################
plot(posterior$samp$k, type="l", main="")


###################################################
### chunk number 45: alt-mixing-nu
###################################################
plot(posterior$samp$nu, type="l", main="")


###################################################
### chunk number 46: alt-mixing-mu
###################################################
plot(posterior$samp$mu, type="l", main="")


###################################################
### chunk number 47: 
###################################################
mean.params <- as.list(apply(posterior$samp, 2, mean))


###################################################
### chunk number 48: 
###################################################
costs <- list(vac=2, death=4, infect=1) 
vacgrid <- list(fracs=seq(0,1.0,0.1), stops=seq(2,init$S0-50,50)) 
alt.optvac <- optvac(init, mean.params, vacgrid, costs, T=time)
alt.best <- getpolicy(alt.optvac)


###################################################
### chunk number 49: alt-optvac
###################################################
plot(alt.optvac)


###################################################
### chunk number 50: 
###################################################
alt.vac.opt <- alt.best[3:4]
alt.MCepi <- MCepi(init, alt.epistep, alt.vac.opt, costs, T=time)


###################################################
### chunk number 51: 
###################################################
getcost(alt.MCepi)


###################################################
### chunk number 52: 
###################################################
alt.MCmanage <- MCmanage(init, alt.epistep, vacgrid, costs, T=time)


###################################################
### chunk number 53: 
###################################################
getcost(alt.MCmanage)


###################################################
### chunk number 54: alt-MCepi-t
###################################################
plot(alt.MCepi, showd=TRUE)


###################################################
### chunk number 55: alt-MCepi-c
###################################################
plot(alt.MCepi, type="costs")


###################################################
### chunk number 56: alt-MCmanage-t
###################################################
plot(alt.MCmanage, showd=TRUE)


###################################################
### chunk number 57: alt-MCmanage-c
###################################################
plot(alt.MCmanage, type="costs")


###################################################
### chunk number 58: 
###################################################
alt.worst <- getpolicy(alt.optvac, which ="worst") 
rbind(alt.best, alt.worst) 


