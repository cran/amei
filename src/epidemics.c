
/* New Code for simulating disease epidemics quickly */

#include "matrix.h"
#include <R.h>
#include <Rmath.h>
#include "epidemics.h"

/* 
   Perform the stochastic forward simulation of the SIR model from Merl & Mangel 2006
   from initial conditions specified by the first entry in the T dimensional S,I,R,D vectors
   containing the numbers of (S)usceptibles, (I)nfectives, (R)ecoveries, (D)eaths.
   
   Note: it may also be reasonable to replace this with the deterministic trajectories of the
   SIRD curves, rather than use stochastic approximations.
*/
void SimulateEpidemic(int* S, int* I, int* R, int* D, int* V, double* C, int T,
		      double b, double k, double nu, double mu,
		      double vacc, int vaccstop,double cvacc, double cdeath, double cinfected,int starttime)
{
  int N0,i,t,itilde,rtilde,dtilde,ztilde,vtilde;
  //int doagain,howmany;
  double *infectionProb,totalRemovalProb,recoveryProb;
  
  N0=S[0]+I[0];
  
  /* precalculate the infections probabilities associated with each number of infecteds */
  infectionProb = new_vector(N0+1);
  for(i=0;i<N0+1;i++)
    infectionProb[i] = 1-exp(k*log(k/(k+b*i)));
  
  /* load the probabilities of removal (meaning either recovery or death), and the relative 
     probabilities of recovery and death within removals */
  totalRemovalProb = 1-exp(-(nu+mu));
  recoveryProb = nu/(nu+mu);
  // mortalityProb = 1- recoveryProb;

  GetRNGstate();

/*   doagain=1; */
/*   howmany=0; */
/*   while(doagain) */
/*     { */
/*       doagain=0; */
  for(t=1;t<T;t++)
    {
      /* 	  if(i==starttime && I[t-1]<=I[0] && howmany<100) */
      /* 	    { */
      /* 	      howmany++; */
      /* 	      if(howmany==100) */
      /* 		{ */
      /* 		  Rprintf("Warning: <1% potential for epidemic\n"); */
      /* 		} */
      /* 	      else */
      /* 		{ */
      /* 		  doagain=1; */
      /* 		  break; */
      /* 		} */
      /* 	    } */
      
      //my_rbinom(1, I[t-1],totalRemovalProb,&ztilde);
      //my_rbinom(1,ztilde,recoveryProb,&rtilde);
      //my_rbinom(1,S[t-1],infectionProb[I[t-1]],&itilde);
      ztilde = rbinom(I[t-1],totalRemovalProb);
      rtilde = rbinom(ztilde,recoveryProb);
      dtilde = ztilde-rtilde;  
      
      // some number are vaccinated, assume that this also curbs infection
      /*       if(S[t-1]>vaccstop) */
      /* 	vtilde=(int) ceil(vacc*(double)S[t=1]); */
      /*       else */
      /* 	vtilde=0; */
      /*       if(S[t-1]>vaccstop) */
      /* 	Rprintf("enacting policy with S=%d\n",S[t-1]); */
      /*       else */
      /* 	Rprintf("not enacting policy with S=%d\n",S[t-1]); */
      
      if(t<starttime)
	{
	  vtilde = 0;
	}
      else
	{
	  vtilde = (S[t-1]>vaccstop) ? (int) ceil(vacc*(double)S[t-1]) : 0;
	}
      
      itilde = rbinom(S[t-1]-vtilde,infectionProb[I[t-1]]);
      
      S[t] = S[t-1] - vtilde - itilde;
      I[t] = I[t-1] + itilde - ztilde;
      R[t] = R[t-1] + rtilde;
      D[t] = D[t-1] + dtilde;
      V[t] = V[t-1] + vtilde;
      C[t] = C[t-1] + cinfected*(double)(I[t-1]-dtilde) + cvacc*(double)vtilde + cdeath*(double)(dtilde);
      
    }
  
  PutRNGstate();	
  free(infectionProb);
}

void RSimulateEpidemic(int* S, int* I, int* R, int* D, int* V, double* C, 
		      int* T, double* b, double* k, double* nu, double* mu,
		       double* vacc, int* vaccstop, double* cvacc, double* cdeath, double* cinfected,int* starttime)
{
  SimulateEpidemic(S,I,R,D,V, C,*T,*b,*k,*nu,*mu,*vacc,*vaccstop,*cvacc,*cdeath,*cinfected,*starttime);
}


void VarStopTimePolicy(int S0,int I0, int T, double b, double k, double nu, double mu, double cvacc, double cdeath, double cinfected,int mcits,double* Vprobs, int nVprobs, int* Vstops,int nVstops, double* EC,int midepidemic,int starttime)
{
  int *S,*I,*R,*D,*V;
  int vaccno,stopno,m,foo,blah;
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
	// if we are considering a vaccination policy that won't be implemented because 
	// the stopnumber is greater than the current number of susceptibles, just use
	// the value calculated for the nonvaccination strategy (vaccno=0).  
	// If we are considering the nonvaccination strategy (vaccno=0), we don't need to
	// repeat simulations for all different stop times, because in each case we're doing nothing.
	if((vaccno>0 && Vstops[stopno]>S0) || (vaccno==0 && stopno>0))
	  {
	    myEC[vaccno][stopno] = myEC[0][0];
	  }
	else
	  {
	    for(m=0;m<mcits;m++)
	      {
		foo=blah=0;
		while(!foo)
		  {
		    SimulateEpidemic(S,I,R,D,V,C,T,b,k,nu,mu,Vprobs[vaccno],Vstops[stopno],cvacc,cdeath,cinfected,starttime);
		    if(!midepidemic)
		      {
			blah++;
			if(blah==100)
			  {
			    foo=1;
			    Rprintf("Warning: <1% chance of an epidemic\n");
			  }
			else
			  {
			    foo = (I[starttime-1]>I[0]);//epicheck(I,T);// || (b<=0. || k<=0.) || Vprobs[vaccno]==1;
			    /* think about this check a bit -- we want to ensure that the epidemic is one that takes off at least
			       a little bit.  we do this by requiring that the size of the epidemic (# infecteds) gets larger than 
			       the initial size at least once during the course of the epidemic.  this could result in an infinite 
			       loop if the epidemic never grows as a result of the vac policy eliminating all susceptibles in one
			       time step.  i think the above takes care of that case */
			  }
		      }
		    else
		      {
			foo=1;
		      }
		  }
		myEC[vaccno][stopno]+=C[T-1];
	      }
	    myEC[vaccno][stopno]/=(double)mcits;
	  }
      }
  free(S);free(I);free(R);free(D);free(V);rm_2ddptrs(myEC);free(C);
  
}

void RVarStopTimePolicy(int* S0,int* I0, int* T, double* b, double* k, double* nu, double* mu, double* cvacc, double* cdeath, double* cinfected,int* mcits,double* Vprobs, int* nVprobs, int* Vstops, int* nVstops, double* EC,int* midepidemic,int* starttime)
{
  VarStopTimePolicy(*S0,*I0,*T,*b,*k,*nu,*mu,*cvacc,*cdeath,*cinfected,*mcits, Vprobs, *nVprobs, Vstops,*nVstops, EC,*midepidemic,*starttime);

}

/*
  check to see if the epidemic actually took off or not, we'll throw away epidemics that don't meet some bare
  requirement for epidemicness
*/
int epicheck(int* I,int T)
{
  int foo,t;
  foo=0;
  for(t=1;t<T;t++)
    if(I[t]>I[0])
      foo++;
  return foo>1;

}
