#include "state.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double inf_prob(State state)
{
  return  1.0 - pow(state.k/(state.k+state.b*(double)state.i), state.k);
}


double postvac_suscept(State state)
{
  return  state.s*(1-((double)state.a)/100);
}


double remaining_inf(State state)
{
  return  state.i*exp(-state.nu-state.mu);
}


double cost_vac(State state)
{
  return state.cv*state.s*(((double)state.a)/100) + state.cd*(((double)state.a)/100)*state.q*state.s;
}


double cost_inf(State state)
{
  return  state.ci*state.i + state.cd*state.i*(1-exp(-state.nu-state.mu))*state.mu/(state.nu+state.mu);
}

