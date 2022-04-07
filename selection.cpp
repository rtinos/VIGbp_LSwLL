/******************************************************************************\
*								 Selection									   *
\******************************************************************************/
#include "defs.h"


/******************************************************************************\
*							Acceptance Criterion: Better		    		   *
\******************************************************************************/
int ACBetter(individual *sol, individual *sol_new)
{
	
	if ( (sol_new->fitness - sol->fitness) > EPS1 ) 
		return (1);
	else
		return (0);

}


/******************************************************************************\
*							Acceptance Criterion: Simulated Annealing  		   *
\******************************************************************************/
int ACSA(individual *sol, individual *sol_new)
{
	double T=0.01;
	
	if ( (sol_new->fitness - sol->fitness) > EPS1 )
		return (1);
	else{
		if ( random_dou() < exp((sol_new->fitness-sol->fitness)/T) )
			return (1);	
		else
			return (0);
	}

}


/******************************************************************************\
*							Acceptance Criterion				    		   *
\******************************************************************************/
void AcceptanceCriterion(individual *sol, individual *sol_new)
{
	int change_flag;
	
	change_flag = ACBetter(sol, sol_new);
	//change_flag = ACSA(sol, sol_new);
	if (change_flag==1){
		for (int gene=0;gene<chrom_size;gene++) 
     		sol->solutionVector[gene] = sol_new->solutionVector[gene];
     	sol->fitness = sol_new->fitness;
	}

}



