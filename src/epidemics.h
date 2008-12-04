#ifndef __EPIDEMICS_H__
#define __EPIDEMICS_H__

void SimulateEpidemic(int* S, int* I, int* R, int* D, int* V, double* C, 
		      int T, double b, double k, double nu, double mu,
		      double vacc, int vaccstop, double cvacc, double cdeath, double cinfected,int starttime);

void RSimulateEpidemic(int* S, int* I, int* R, int* D, int* V, double* C, 
		      int* T, double* b, double* k, double* nu, double* mu,
		       double* vacc, int* vaccstop, double* cvacc, double* cdeath, double* cinfected,int* starttime);

void VarStopTimePolicy(int S0,int I0, int T, double b, double k, double nu, double mu, double cvacc, double cdeath, double cinfected,int mcits,double* Vprobs, int nVprobs, int* Vstops,int nVstops, double* EC,int midepidemic,int starttime);

void RVarStopTimePolicy(int* S0,int* I0, int* T, double* b, double* k, double* nu, double* mu, double* cvacc, double* cdeath, double* cinfected,int* mcits,double* Vprobs, int* nVprobs, int* Vstops, int* nVstops, double* EC,int* midepidemic, int* starttime);

int epicheck(int* I,int T);

#endif
