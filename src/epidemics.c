
/* C Code for simulating disease epidemics quickly */
#include "matrix.h"
#include <R.h>
#include <Rmath.h>
#include "epidemics.h"

/* 
   Perform the stochastic forward simulation of the SIR model from
   Merl & Mangel 2006 from initial conditions specified by the first
   entry in the T dimensional S,I,R,D vectors containing the numbers
   of (S)usceptibles, (I)nfectives, (R)ecoveries, (D)eaths.
   
   Note: it may also be reasonable to replace this with the
   deterministic trajectories of the SIRD curves, rather than use
   stochastic approximations.

   Parameters are:
   S: # susceptibles at each time step
   I: # infecteds at each time step
   R: # recovereds at each time step
   D: # deaths at each time step
   V: # vaccinated at each time step
   C: acumulated cost at each time step
   T: number time steps to simulate
   b: transmission rate
   k: dispersion parameter
   nu: recovery rate
   mu: mortality rate
   vacc: percentage to vaccinate
   vaccstop: number of susceptibles below which to stop vaccination
   cvacc: cost per vaccination
   cdeath: cost per death
   cinfected: cost per infected per time unit
   starttime: time to begin vaccinations
*/
void SimulateEpidemic(int* S, int* I, int* R, int* D, int* V, double* C, int T,
		      double b, double k, double nu, double mu,
		      double vacc, int vaccstop,double cvacc, double cdeath, 
		      double cinfected,int starttime)
{
  int N0,t,itilde,rtilde,dtilde,ztilde,vtilde;
  double totalRemovalProb,recoveryProb;
  
  N0=S[0]+I[0];
  
  /* precalculate the infections probabilities associated with each
     number of infecteds */
  /* infectionProb = new_vector(N0+1);
     for(i=0;i<N0+1;i++)
     infectionProb[i] = 1-exp(k*log(k/(k+b*i))); */
  
  /* load the probabilities of removal (meaning either recovery or
     death), and the relative probabilities of recovery and death
     within removals */
  totalRemovalProb = 1-exp(-(nu+mu));
  recoveryProb = nu/(nu+mu);

  GetRNGstate();

  for(t=1;t<T;t++)
    {
      ztilde = rbinom(I[t-1],totalRemovalProb);
      rtilde = rbinom(ztilde,recoveryProb);
      dtilde = ztilde-rtilde;  
            
      if(t<starttime)
	{
	  vtilde = 0;
	}
      else
	{
	  vtilde = (S[t-1]>vaccstop) ? (int) ceil(vacc*(double)S[t-1]) : 0;
	}
      
      /* Thanks to Michael Hohle */
      /* itilde = rbinom(S[t-1]-vtilde,infectionProb[I[t-1]]); */
      itilde = rbinom(S[t-1]-vtilde, 1-exp(k*log(k/(k+b*I[t-1]))));
      
      S[t] = S[t-1] - vtilde - itilde;
      I[t] = I[t-1] + itilde - ztilde;
      R[t] = R[t-1] + rtilde;
      D[t] = D[t-1] + dtilde;
      V[t] = V[t-1] + vtilde;
      C[t] = C[t-1] + cinfected*(double)(I[t-1]-dtilde) + cvacc*(double)vtilde + 
	cdeath*(double)(dtilde);
      
    }
  
  PutRNGstate();	
  /* free(infectionProb); */
}

// R interface to SimulateEpidemics()
void RSimulateEpidemic(int* S, int* I, int* R, int* D, int* V, double* C, 
		      int* T, double* b, double* k, double* nu, double* mu,
		       double* vacc, int* vaccstop, double* cvacc, 
		       double* cdeath, double* cinfected,int* starttime)
{
  SimulateEpidemic(S,I,R,D,V, C,*T,*b,*k,*nu,*mu,*vacc,*vaccstop,
		   *cvacc,*cdeath,*cinfected,*starttime);
}

/*
  Calculates the variable stop time policy for initial conditions: 
  S0: initial number of susceptibles 
  I0: initial number of infecteds
  T: number of time steps to consider 
  b: transmission rate 
  k: dispersion parameter 
  nu: recovery rate 
  mu: mortality rate 
  cvacc: cost per vaccination 
  cdeath: cost per death 
  cinfected: cost per infected per time unit 
  mcits: number mc iterations to average over when calculating expected cost
  vprobs: vector of probabilities representing possible vaccination policies 
  nVprobs: length of vprobs vector 
  vstops: vector of thresholds representing possible stop times
  nVstops: length of vstops vector
  EC: matrix of expected costs for each vprob x vstop pair 
  midepidemic: boolean indicating whether this is a policy being calculated
    to implement in the middle of an epidemic
 starttime: time step during which to begin vaccination
*/
void VarStopTimePolicy(int S0,int I0, int T, double b, double k, double nu, 
		       double mu, double cvacc, double cdeath, double cinfected,
		       int mcits,double* Vprobs, int nVprobs, int* Vstops,int nVstops, 
		       double* EC,int midepidemic,int starttime)
{
  int *S,*I,*R,*D,*V;
  int vaccno,stopno,m,nvalid,ntries;
  double **myEC, *C;

  S = new_ivector(T); // susceptibles
  S[0] = S0;
  I = new_ivector(T); // infecteds
  I[0] = I0;
  R = new_ivector(T); // recovereds
  R[0] = 0;
  D = new_ivector(T); // deaths
  D[0] = 0;
  V = new_ivector(T); // vaccinated
  V[0] = 0;
  C = new_vector(T);  // costs incurred
  C[0] = 0;

  zerov(EC,nVprobs*nVstops);  // expected costs per policy
  myEC = new_2ddptrs(nVprobs,nVstops,EC);
  
  for(vaccno=0;vaccno<nVprobs;vaccno++)
    for(stopno=0;stopno<nVstops;stopno++)
      {
	/*
	if we are considering a vaccination policy that won't be
	implemented because the stopnumber is greater than the current
	number of susceptibles, just use the value calculated for the
	nonvaccination strategy (vaccno=0).  If we are considering the
	nonvaccination strategy (vaccno=0), we don't need to repeat
	simulations for all different stop times, because in each case
	we're doing nothing.
	*/
	if((vaccno>0 && Vstops[stopno]>S0) || (vaccno==0 && stopno>0))
	  {
	    myEC[vaccno][stopno] = myEC[0][0];
	  }
	else
	  {
	    for(m=0;m<mcits;m++)
	      {
		nvalid=ntries=0;
		while(!nvalid)
		  {
		    SimulateEpidemic(S,I,R,D,V,C,T,b,k,nu,mu,Vprobs[vaccno],
				     Vstops[stopno],cvacc,cdeath,cinfected,starttime);
		    if(!midepidemic)
		      {
			ntries++;
			if(ntries==100)
			  {
			    nvalid=1;
			    Rprintf("Warning: <1% chance of an epidemic\n");
			  }
			else
			  {
			    nvalid = (I[starttime-1]>I[0]);
			    
			  }
		      }
		    else
		      {
			nvalid=1;
		      }
		  }
		myEC[vaccno][stopno]+=C[T-1];
	      }
	    myEC[vaccno][stopno]/=(double)mcits;
	  }
      }
  free(S);free(I);free(R);free(D);free(V);rm_2ddptrs(myEC);free(C);
  
}

// R interface to VarStopTimePolicy()
void RVarStopTimePolicy(int* S0,int* I0, int* T, double* b, double* k, double* nu, 
			double* mu, double* cvacc, double* cdeath, double* cinfected,
			int* mcits,double* Vprobs, int* nVprobs, int* Vstops, int* nVstops,
			double* EC,int* midepidemic,int* starttime)
{
  VarStopTimePolicy(*S0,*I0,*T,*b,*k,*nu,*mu,*cvacc,*cdeath,*cinfected,*mcits,
		    Vprobs, *nVprobs, Vstops,*nVstops, EC,*midepidemic,*starttime);

}

/*
  check to see if the epidemic actually took off or not, we'll throw
  away epidemics that don't meet some bare requirement for
  epidemicness
*/
int epicheck(int* I,int T)
{
  int nvalid,t;
  nvalid=0;
  for(t=1;t<T;t++)
    if(I[t]>I[0])
      nvalid++;
  return nvalid>1;

}
