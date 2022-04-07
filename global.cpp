#include "defs.h"

// Global variables and parameters
int chrom_size;										// size of the solutionVector (dimension of the problem) 
int gen;											// generation
int n_instances = 10;								// number of instances (for NK landscape problem)
int n_runs_per_instance = 10;						// number of runs for each instance (for NK landscape problem)
int save_extra_info_flag = 1;						// flag for saving extra information
int perturbation_type;								// Perturbation method: 1. kick solution (SRP); 2. estimated VIG based (black box VIGbp); 3. real VIG based (gray box VIGbp); 4. same as 1, but random size (RPN); 5. same as 1, but using rate (RPR) 
int total_edges_VIG=0;								// Total number of edges in the VIG (used for statistics in perturbation 2)
long int max_gen;									// maximum number of generations (when stop criterion is time, it is used for saving some statistics data)
// Vectors
int *vsort_aux;										// auxiliar sorted vector of integers (used in different methods)
int *file_gen;										// data to be stored: number of generations								
double *file_best_fitness, *time_run;				// data to be stored: best fitness, runtime
double *file_best_fitness_gen;						// data to be stored: best fitness over the generations 
double *file_n_edges_eVIG_gen;						// data to be stored: number of edges of the eVIG
double *file_n_iter_ls;								// data to be stored: number of iterations of Local Search
double *file_n_best_fitness_improv;					// data to be stored: number of improvements of the best solution
double *file_n_escapes;								// data to be stored: number of escapes from local optimum
double *file_HD_LO;									// data to be stored: Hamming distance between two local optima
double *file_HD_pert;								// data to be stored: Hamming distance between solution before and after perturbation
double *file_fl_pert;								// data to be stored: number of subfunctions changed by perturbation
// Matrices	
int **File_best_ind;								// data to be stored: best individual

