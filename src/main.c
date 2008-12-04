/****************************************************************/
/* FILE: main.c		 				        */
/*                                                              */
/* AUTHOR: Leah Johnson (from R code by Dan Merl)	       	*/
/* DATE:   October 11, 2004	       				*/
/*     								*/
/* DESCRIPTION: This program contains a function which uses     */
/*		policy iteration to solve the optimal           */
/*              vaccination policy problem by dynamic		*/
/*              programming.                                    */
/*	       	The scenario here is that of an infectious disease  */
/*    		unfolding within a population of initial size N.    */
/*		For each division of the N-size population into S   */
/*		initial susceptibles and I intitial infecteds, the  */
/*    		function determines the expected number of deaths   */
/*    		under optimal vaccination, and the fraction of the  */
/*    		population to vaccinate.			 */
/*								 */
/* CONSTANTS:	mcits   100  - number of draws to use for monte  */
/*			       carlo estimation of expectation   */
/*			       over itilde			 */
/*		thresh  .01  - difference threshold used to	 */
/*			       measure convergence.  When the 	 */
/*			       sum of the squares of the	 */
/*			       differences of the elements in	 */
/*			       Vnew and Vold is less than thresh,*/
/*			       convergence is assumed.		 */
/*		maxiters 10  - max number of iterations allowed, */
/*			       regardless of convergence. Making */
/*			       this a small number can prevent   */
/*			       convergence, but will also keep   */
/*			       the code from running forever.	 */
/*			 					 */
/* VARIABLES:   The underlying SI model contains the following	 */
/*		parameters, saved into State:			 */
/*								 */
/*		float   N	- initial population size	 */
/*			b	- transmission parameter	 */
/*			k	- interaction parameter		 */
/*			m	- mortality rate		 */
/*			n	- recovery rate			 */
/*			q	- prob of death by vaccination	 */
/*			c	- cost of a single infection	 */
/*			cv	- cost of a single vaccination	 */
/*			cd	- "cost" of single death         */
/*								 */
/*		Other variables:				 */
/*								 */
/*		int	s, i, j, a  - counters (in State, except j)     */
/*		double  *alphavals  - stores the expected number of     */
/*		      		      deaths under each vaccination	*/
/*			       	      policy. In this case we assume    */
/*				      101 possible policies:	        */
/*				      0, 0.01, 0.02, ... , .99, 1	*/
/*				      (fraction of suscept vaccinated)  */
/*		double  *itilde - sample from the distribution of       */
/*			          new infecteds. (here binomial)        */
/*		int	*index1	- indices calculated for the mc	        */
/*			*index2	  for Vold or Vnew		        */
/*		double  *sum	- sum of entries from Vold or Vnew      */
/*			prob	- prob of inf (from inf_prob(state))    */
/*			suscept - # suscept after vaccination (from     */
/*				  postvac_suscept(state))		*/
/*									*/
/*		Execution creates the following variables:		*/
/*									*/
/*		double   **Vnew - matrix containing expected number     */
/*				  of deaths under optimal		*/
/*				  vaccination. Rows correspond to       */
/*				  infecteds, columns to susceptibles    */
/*		double   **Vold - matrix containing previous		*/
/*				  iteration of Vnew			*/
/*		double   **alpha -matrix containing optimal fraction    */ 
/*				  of susceptibles to vaccinate. Rows    */
/*				  correspond to infecteds, columns      */
/*				  to susceptibles			*/
/*		double   **alpha_num -matrix containing optimal		*/ 
/*				  number of suscept to vaccinate.       */
/*				  rows -> infecteds, col - >suscept     */
/*		double   iters  - number of iterations required for     */
/*				  the algorithm to converge		*/
/*									*/
/*									*/
/* USAGE:								*/
/************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "ssdiff.h"
//#include "rand_draws.h"
#include "state.h"
#include <R.h>
#include <Rmath.h>

//#define PI		3.14159265358979323846
#define mcits		100
#define thresh		10e-3
#define maxiters	100

void newvacpolicy(int *N_in, double *b_in, double* k_in, double *mu_in, 
		  double *nu_in, double *q_in, double *ci_in, double *cv_in, 
		  double *cd_in, double *thresh_in, int *mcits_in, 
		  int *maxiters_in,
		  double *V_in, double *V_out, double *alpha_out, 
		  double *alpha_num_out, int* numiters);

double **initializeVold(State state, int *index1, int *index2, double *itilde,
			int mymcits);

/* int main (void /\*int argc, const char * argv[]*\/) { */
/*   float	b, k, mu, nu, q, ci, cv, cd, mythresh; */
/*   int	N, mymaxiters; */
/*   int   mymcits[2]; */
  
/*   printf("Please input the following parameter values:\n\n"); */
/*   printf("\nInitial population size, N: "); */
/*   scanf("%d", &N); */
/*   printf("\nTransmission parameter, b: "); */
/*   scanf("%f", &b); */
/*   printf("\nInteraction parameter, k: "); */
/*   scanf("%f", &k); */
/*   printf("\nMortality rate, mu: "); */
/*   scanf("%f", &mu); */
/*   printf("\nRecovery rate, nu: "); */
/*   scanf("%f", &nu); */
/*   printf("\nProbibility of death from vaccination, q: "); */
/*   scanf("%f", &q); */
/*   printf("\nCost of infection, in dollars, ci: "); */
/*   scanf("%f", &ci); */
/*   printf("\nCost of vaccination, in dollars, cv: "); */
/*   scanf("%f", &cv); */
/*   printf("\nCost of a single death, in dollars, cd: "); */
/*   scanf("%f", &cd); */
  
/*   printf("\nThese are the parameters you entered: \n" \ */
/* 	 "N = %d \n b = %f \n k = %f \n mu = %f \n" \ */
/* 	 "nu = %f \n q= %f ci = %f cv = %f cd = %f\n\n",  */
/* 	 N, b, k, mu, nu, q, ci, cv, cd); */
  
/*   mythresh = thresh; */
/*   mymcits[0]= mcits; */
/*   mymcits[1]= 10*10*10*mcits; */
/*   mymaxiters = maxiters; */

/*   newvacpolicy(&N, (double*) &b, (double*) &k, (double*) &mu,  */
/* 	       (double*) &nu, (double*) &q,(double*) &ci,  */
/* 	       (double*) &cv, (double*) &cd, (double*) &mythresh,  */
/* 	       (int*) &mymcits, (int*) &mymaxiters,  */
/* 	       NULL, NULL, NULL, NULL, NULL); */
  
/*   return 0; */
/* } */


void newvacpolicy(int *N_in, double *b_in, double* k_in, double *mu_in, 
		  double *nu_in, double *q_in, double *ci_in, double *cv_in, 
		  double *cd_in, double *thresh_in, int *mcits_in, 
		  int *maxiters_in, double *V_in, 
		  
		  double *V_out, double *alpha_out, 
		  double *alpha_num_out, int *numiters)
{
  double	sum, diff, prob, lastdiff;
  int		j, iters, suscept;
  unsigned int which, M;
  double	**Vold, **Vnew, **alpha, **alpha_num;
  double	*alphavals, *itilde;
  int		*index1, *index2;

  lastdiff= 1e300*1e300;
	
  State state;
  state.N = *N_in; state.b = *b_in; state.k = *k_in; state.mu = *mu_in; 
  state.nu = *nu_in, state.q = *q_in, state.ci=*ci_in, state.cv=*cv_in; 
  state.cd=*cd_in;	
  
  Vnew = new_zero_matrix(state.N+1,state.N+1);
  alphavals = new_zero_vector(101);
  alpha = new_zero_matrix(state.N+1, state.N+1);
  alpha_num = new_zero_matrix(state.N+1, state.N+1);
  
  index1 = (int*) malloc(sizeof(int)* (mcits_in[1]));
  index2 = (int*) malloc(sizeof(int)* (mcits_in[1]));
  M = (unsigned int)(state.N+1.0);
  itilde = new_vector(mcits_in[1]);
  
  if( V_in == NULL){
    Vold = initializeVold(state, index1, index2, itilde, mcits_in[0]);
  }else{
    Vold = new_matrix(state.N+1,state.N+1);
    dupv( Vold[0], V_in , (state.N+1)*(state.N+1) ); 
  }
	
  iters = 0;
  
  while( 1 ){
    for( state.s = 0; state.s <= (state.N-1); state.s++ ){
      for( state.i = 1; state.i <= (state.N - state.s); state.i++ ){
	for( state.a = 0; state.a <= 100; state.a++ ){
	  
	  /* sample from the distribution of new infecteds */
	  prob = inf_prob(state); 
	  suscept = floor(postvac_suscept(state)); 
	  //for( j = 0; j < mcits_in[0]; j++ ) itilde[j] = bnldev( prob, suscept );
	  for( j = 0; j < mcits_in[0]; j++ ) itilde[j] = rbinom(suscept,prob);

	  /* store the expected number of deaths under for that vaccination  */
	  /* policy using Monte Carlo simulation to compute the expectation  */
	  /* of Vold over itilde.	  			             */
	  /* mean(Vold[floor(i*exp(-nu-mu)+itilde+1),                        */
	  /* floor(s*(1-a/100)-itilde+1)]))                                  */
	  for( j = 0; j < mcits_in[0]; j++ ){
	    index1[j] = (int) (remaining_inf(state)+itilde[j]);
	    index2[j] = (int) (postvac_suscept(state)-itilde[j]);
	  }
	  
	  sum = 0;
	  for( j = 0; j < mcits_in[0]; j++)sum = sum + Vold[index1[j]][index2[j]];
	  
	  /* alphavals[a+1]<-(a/100*q*s+i*(1-exp(-nu-mu))*mu/(nu+mu)+ 	    */
	  /* mean(Vold[floor(i*exp(-nu-mu)+itilde+1),                       */
	  /* floor(s*(1-a/100)-itilde+1)]))                                 */
	  alphavals[state.a] = cost_vac(state) + cost_inf(state) + sum/(mcits_in[0]);
	}
	
	/* save the optimal (least) expected number of deaths or least cost */
	Vnew[state.i][state.s] = min(alphavals, 101, &which);
	
	/* save the corresponding vaccination policy */
	alpha[state.i][state.s]=((double) which)/100;
	
	/* save the vaccination policy in terms of # of people to vaccinate */
	alpha_num[state.i][state.s]=state.s*alpha[state.i][state.s];
	
      }
    }
    
    
    /* if we haven't reached convergence or gone too many iterations,   */
    /* increase the iteration count, and copy Vnew into Vold.  	     	*/
    
    diff = ssdiff( Vnew, Vold, M, M );
    printf("iter=%d, diff = %g, thresh = %g, mcits = %d \n", 
	   iters, diff, *thresh_in, mcits_in[0]);
    
    if( diff > (*thresh_in) && iters < *maxiters_in ){  

      if( diff > lastdiff ){
	mcits_in[0]*=10;
	if( mcits_in[0] > mcits_in[1]) break;
	printf("oscillating, incremented mcits, %d \n", mcits_in[0]);
	lastdiff = 1e300*1e300;
      }else{
	lastdiff = diff;
      }

      iters = iters+1;
      swap_matrix(Vold, Vnew, M, M);
      
      continue;
    }
    /* if we've reached convergence or gone too many iterations, */
    /*   break out of while loop */
    else{
      if(mcits_in[0]<= 100){
	mcits_in[0]*=10;
	if( mcits_in[0] > mcits_in[1]) break;
	printf("premature convergance, incremented mcits, %d \n", mcits_in[0]);
	lastdiff = 1e300*1e300;
	iters = iters+1;
	swap_matrix(Vold, Vnew, M, M);
	continue;
      }
      break;
    }
  }
  
  /* at this point, Vnew should have the optimal expected deaths and   	*/
  /* alpha should have the optimal vaccination fractions       		*/
  
  /*printf("\n number of iterations: %d \n\n", iters);*/
  
  /*printf("Vnew = \n");
    printMatrix(Vnew, N+1, N+1, stdout);
    printf("\nalpha = \n");
    printMatrix(alpha, N+1, N+1, stdout);*/

  *numiters = iters;
  
  if(V_out) dupv(V_out, Vnew[0], (state.N+1)*(state.N+1));
  if(alpha_out) dupv(alpha_out, alpha[0], (state.N+1)*(state.N+1));
  if(alpha_num_out) dupv(alpha_num_out, alpha_num[0], (state.N+1)*(state.N+1));
  
  delete_matrix(alpha);
  delete_matrix(alpha_num);
  delete_matrix(Vold);
  delete_matrix(Vnew);
  free(index1);
  free(index2);
  free(itilde);
  free(alphavals);
}



double **initializeVold(State state, int *index1, int *index2, double *itilde,
			int mymcits)
{
  double        **Vold;
  double	sum, prob;
  int		j;

  Vold = new_zero_matrix(state.N+1,state.N+1);
  
  for( state.s = 0; state.s <= (state.N-1); state.s++ ){
    for( state.i = 1; state.i <= (state.N - state.s); state.i++ ){   
      /* sample from the distribution of new infecteds*/
      
      prob = inf_prob(state); 
      //      for( j = 0; j < mymcits; j++ ) itilde[j] = bnldev(prob, state.s);
      for( j = 0; j < mymcits; j++ ) itilde[j] = rbinom(state.s, prob);

      /* store the expected number of deaths without vaccination using  */
      /* Monte Carlo simulation to compute the expectation of Vold over */
      /* itilde.							*/
      /* mean(Vold[floor(i*exp(-nu-mu)+itilde+1), floor(s-itilde+1)])	*/
      for( j = 0; j < mymcits; j++ ){
	index1[j] = (int) (remaining_inf(state) + itilde[j]);
	index2[j] = (int) (state.s-itilde[j]);
      }
      
      sum = 0;
      for( j = 0; j < mymcits; j++) sum = sum + Vold[index1[j]][index2[j]];
      
      /* recursively compute the expected number of deaths using Monte  */ 
      /* Carlo estimation for the expectation of Vold over itilde	*/
      
      Vold[state.i][state.s] = cost_inf(state) + sum/(mymcits);
    }
  }
  return Vold;
}
