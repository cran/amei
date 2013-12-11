### R code from vignette source 'amei.Rnw'

###################################################
### code chunk number 1: amei.Rnw:21-23
###################################################
library(amei)
options(width=70, prompt="R> ", digits=5)


###################################################
### code chunk number 2: pnas.iRnw:3-4
###################################################
set.seed(12345)


###################################################
### code chunk number 3: pnas.iRnw:15-18
###################################################
tp <- list(b=0.00218, k=10, nu=0.4, mu=0) 
init <- list(S0=762, I0=1, R0=0, D0=0) 
costs <- list(vac=2, death=4, infect=1)


###################################################
### code chunk number 4: pnas.iRnw:34-35
###################################################
vac <- list(frac=0, stop=0)


###################################################
### code chunk number 5: pnas.iRnw:49-50
###################################################
init.MCepi <- MCepi(init, tp, vac, costs)


###################################################
### code chunk number 6: epis
###################################################
plot(init.MCepi)


###################################################
### code chunk number 7: costs
###################################################
plot(init.MCepi, type="costs")


###################################################
### code chunk number 8: pnas.iRnw:83-84
###################################################
vacgrid <- list(fracs=seq(0,1.0,0.1),stops=seq(2,init$S0-75,75)) 


###################################################
### code chunk number 9: pnas.iRnw:92-93
###################################################
out.optvac <-optvac(init, tp, vacgrid, costs) 


###################################################
### code chunk number 10: pnas.iRnw:98-101
###################################################
best <-getpolicy(out.optvac) 
worst <-getpolicy(out.optvac, which="worst")
rbind(best, worst)


###################################################
### code chunk number 11: optvac
###################################################
plot(out.optvac) 


###################################################
### code chunk number 12: pnas.iRnw:131-133
###################################################
vac.opt <- best[3:4]
opt.MCepi <- MCepi(init, tp, vac.opt, costs)


###################################################
### code chunk number 13: episov
###################################################
plot(opt.MCepi)


###################################################
### code chunk number 14: costsov
###################################################
plot(opt.MCepi, type="costs")


###################################################
### code chunk number 15: pnas.iRnw:163-164
###################################################
getvac(opt.MCepi)


###################################################
### code chunk number 16: pnas.iRnw:176-180
###################################################
T <- length(opt.MCepi$Median$C)
optC <- getcost(opt.MCepi)
initC <- getcost(init.MCepi)
data.frame(rbind(initC,optC), row.names=c("init", "opt"))


###################################################
### code chunk number 17: pnas.iRnw:188-189
###################################################
opt.MCepi


###################################################
### code chunk number 18: pnas.iRnw:235-236
###################################################
out.man <- manage(init, epistep, vacgrid, costs) 


###################################################
### code chunk number 19: epi
###################################################
plot(out.man)


###################################################
### code chunk number 20: cost
###################################################
plot(out.man, type="cost")


###################################################
### code chunk number 21: pnas.iRnw:263-264
###################################################
getcost(out.man)


###################################################
### code chunk number 22: params
###################################################
true <- as.list(formals(epistep)$true)
plot(out.man, type="params",tp=tp) 


###################################################
### code chunk number 23: pnas.iRnw:290-291
###################################################
out.man


###################################################
### code chunk number 24: pnas_mcmanage.iRnw:1-2
###################################################
set.seed(12345)


###################################################
### code chunk number 25: pnas_mcmanage.iRnw:16-17
###################################################
out.MCmanage <- MCmanage(init, epistep, vacgrid, costs)


###################################################
### code chunk number 26: MCmanepis
###################################################
plot(out.MCmanage)


###################################################
### code chunk number 27: MCmancosts
###################################################
plot(out.MCmanage, type="costs")


###################################################
### code chunk number 28: pnas_mcmanage.iRnw:46-47
###################################################
getvac(out.MCmanage)


###################################################
### code chunk number 29: MCmanfracs
###################################################
plot(out.MCmanage, type="fracs")


###################################################
### code chunk number 30: MCmanstops
###################################################
plot(out.MCmanage, type="stops")


###################################################
### code chunk number 31: pnas_mcmanage.iRnw:71-76
###################################################
cinit <- getcost(init.MCepi)
copt <- getcost(opt.MCepi)
cman <- getcost(out.MCmanage)
data.frame(rbind(cinit, copt, cman), 
           row.names=c("init", "opt", "man"))


###################################################
### code chunk number 32: pnas_mcmanage.iRnw:91-92
###################################################
bad <- list(b=0.001, k=10, nu=0.9, mu=0)


###################################################
### code chunk number 33: pnas_mcmanage.iRnw:96-99
###################################################
costs.bad <- optvac(init, bad, vacgrid, costs)
pol.bad <- getpolicy(costs.bad)
pol.bad


###################################################
### code chunk number 34: pnas_mcmanage.iRnw:104-107
###################################################
bad.MCepi <- MCepi(init, tp, pol.bad[3:4], costs)
cbad <- getcost(bad.MCepi)
cbad


###################################################
### code chunk number 35: alt_trans.iRnw:6-7
###################################################
set.seed(12345)


###################################################
### code chunk number 36: alt_trans.iRnw:67-96
###################################################
alt.epistep <- 
function(SIR, last=list(rem=0, rec=0, infect=0, dead=0, Z=0),
         tp=list(a = 0.05, mu = 0.05, nu = 0.1, m = 0.4, 
	 rho = 200, C = 500))
{
  ## calculate the infection probability based on the
  ## reservoir, and randomly infect susceptibles
  Z <- last$Z
  fz <- Z/(Z+tp$C)
  pi <- 1 - exp(-tp$a * fz)
  infect <- rbinom(1, SIR$S, pi)

  ## update recovereds and deaths
  pr <- 1 - exp(-tp$nu)
  rec <- rbinom(1,SIR$I,pr)
  pd <- 1 - exp(-tp$mu)
  dead <- rbinom(1, SIR$I-rec, pd)

  ## reservoir dynamics
  pz <- 1 - exp(-tp$m)
  dz <- rbinom(1, Z, pz)
  bz <- round(SIR$I*tp$rho)
  Z <- Z - dz + bz

  ## the returned list is passed in as "last" in a
  ## subsequent call to this "epistep" function
  return(list(rem=(rec+dead), rec=rec, infect=infect, 
              dead=dead, Z=Z))
}


###################################################
### code chunk number 37: alt_trans.iRnw:108-113
###################################################
init1 <- list(S0=150, I0=1, R0=0, D0=0)
tp <- list(a=0.1, mu=0.0, nu=0.3, m=50, rho=500, C=500)
alt.epistep1 <- alt.epistep
formals(alt.epistep1)$tp <- tp
out.alt<- manage(init1, alt.epistep1, NULL, NULL, T=80)


###################################################
### code chunk number 38: alt_trans.iRnw:123-127
###################################################
tp <- list(a=0.1, mu=0.0, nu=0.3, m=0.001, rho=500, C=500)
alt.epistep2 <- alt.epistep
formals(alt.epistep2)$tp <- tp
out.alt2 <- manage(init1, alt.epistep2, NULL, NULL, T=80)


###################################################
### code chunk number 39: alt
###################################################
plot(out.alt,showv=FALSE)


###################################################
### code chunk number 40: alt2
###################################################
plot(out.alt2,showv=FALSE)


###################################################
### code chunk number 41: alt-params
###################################################
plot(out.alt, type="params", showd=TRUE)


###################################################
### code chunk number 42: alt2-params
###################################################
plot(out.alt2, type="params", showd=TRUE)


###################################################
### code chunk number 43: alt_trans.iRnw:214-217
###################################################
init <- list(S0=600, I0=1, R0=0, D0=0)
time <- 80
posterior <- manage(init, alt.epistep, NULL, NULL, T=time, bkrate=100)


###################################################
### code chunk number 44: alt-mixing-b
###################################################
plot(log(posterior$samp$b), type="l", main="",ylab=expression(b))


###################################################
### code chunk number 45: alt-mixing-k
###################################################
plot(posterior$samp$k, type="l", main="",ylab=expression(k))


###################################################
### code chunk number 46: alt-mixing-nu
###################################################
plot(posterior$samp$nu, type="l", main="",ylab=expression(nu))


###################################################
### code chunk number 47: alt-mixing-mu
###################################################
plot(posterior$samp$mu, type="l", main="",ylab=expression(mu))


###################################################
### code chunk number 48: alt_trans.iRnw:246-247
###################################################
mean.params <- as.list(apply(posterior$samp, 2, mean))


###################################################
### code chunk number 49: alt_trans.iRnw:252-256
###################################################
costs <- list(vac=2.5, death=4, infect=1) 
vacgrid <- list(fracs=seq(0,1.0,0.1), stops=seq(2,init$S0-50,50)) 
alt.optvac <- optvac(init, mean.params, vacgrid, costs, T=time)
alt.best <- getpolicy(alt.optvac)


###################################################
### code chunk number 50: alt-optvac
###################################################
plot(alt.optvac)


###################################################
### code chunk number 51: alt_trans.iRnw:274-276
###################################################
alt.vac.opt <- alt.best[3:4]
alt.MCepi <- MCepi(init, alt.epistep, alt.vac.opt, costs, T=time)


###################################################
### code chunk number 52: alt_trans.iRnw:278-279
###################################################
getcost(alt.MCepi)


###################################################
### code chunk number 53: alt_trans.iRnw:286-287
###################################################
alt.MCmanage <- MCmanage(init, alt.epistep, vacgrid, costs, T=time)


###################################################
### code chunk number 54: alt_trans.iRnw:289-290
###################################################
getcost(alt.MCmanage)


###################################################
### code chunk number 55: alt-MCepi-t
###################################################
plot(alt.MCepi, showd=TRUE)


###################################################
### code chunk number 56: alt-MCepi-c
###################################################
plot(alt.MCepi, type="costs")


###################################################
### code chunk number 57: alt-MCmanage-t
###################################################
plot(alt.MCmanage, showd=TRUE)


###################################################
### code chunk number 58: alt-MCmanage-c
###################################################
plot(alt.MCmanage, type="costs")


###################################################
### code chunk number 59: alt_trans.iRnw:343-345
###################################################
alt.worst <- getpolicy(alt.optvac, which ="worst") 
rbind(alt.best, alt.worst) 


