#ifndef __STATE_H__
#define __STATE_H__ 

/*typedef enum METHOD {DOLS=101, DEATHS=102, DD=103} METHOD;*/

typedef struct state{
			int		a;
			int		s;
			int		i;
			int		N;
			double  b;
			double  k;
			double  mu;
			double  nu;
			double  q;
			double  ci;
			double  cv;
			double  cd;
		} State;
	
double inf_prob(State state);
double postvac_suscept(State state);
double remaining_inf(State state);
double cost_vac(State state);
double cost_inf(State state);

#endif
