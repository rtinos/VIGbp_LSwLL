/******************************************************************************\
*								 Statistics						 			   *
\******************************************************************************/
#include "defs.h"


/******************************************************************************\
*				  Record Solutions											   *
\******************************************************************************/
void statistics(individual *sol, int n_run, long int n_iter_ls)
{	

	// Save best solution of the run and information about improvement and number of iterations of LS
	if (gen==0 || (sol->fitness - file_best_fitness[n_run])>EPS1){
		file_best_fitness[n_run] = sol->fitness;
		for (int gene=0;gene<chrom_size;gene++) 
     		File_best_ind[n_run][gene]=sol->solutionVector[gene];
     	if (gen!=0 && gen<max_gen)	
			file_n_best_fitness_improv[n_run]+=1.0;		// compute only in max_gen generations
	}
	file_n_iter_ls[n_run]+=n_iter_ls;
	
	// Save best fitness along the generations: only for the first run
	if (save_extra_info_flag==1){			
		if (gen<max_gen)			
			file_best_fitness_gen[gen]+=file_best_fitness[n_run];				
	}
			
}


