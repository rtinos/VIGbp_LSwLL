/******************************************************************************\
*								 Definitions							 *
\******************************************************************************/
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <list> 
#define CHAR_LEN 1000
#define EPS1 1.0e-12
#define EPS2 1.0e-8

using namespace std; 

// Data structures
typedef int variableType; 									// type variableType
typedef struct {
			variableType *solutionVector;								
			double fitness;									
} individual;												// data structure individual


// Global variables and parameters
extern int chrom_size;										// size of the solutionVector (dimension of the problem) 
extern int gen;												// generation (iteration)
extern int n_instances;										// number of instances (for NK landscape problem)
extern int n_runs_per_instance;								// number of runs for each instance (for NK landscape problem)
extern int save_extra_info_flag;							// flag for saving extra information
extern int perturbation_type;								// Perturbation method: 1. kick solution (SRP); 2. estimated VIG based (black box VIGbp); 3. real VIG based (gray box VIGbp); 4. same as 1, but random size (RPN); 5. same as 1, but using rate (RPN) 
extern int total_edges_VIG;									// Total number of edges in the VIG (used for statistics in perturbation 2)exter
extern long int max_gen;									// maximum number of generations (when stop criterion is time, it is used for saving some statistics data)
// Vectors
extern int *vsort_aux;										// auxiliar sorted vector of integers (used in different methods)
extern int *file_gen;										// data to be stored: number of generations								
extern double *file_best_fitness, *time_run;				// data to be stored: best fitness, runtime
extern double *file_best_fitness_gen;						// data to be stored: best fitness over the generations 
extern double *file_n_edges_eVIG_gen;						// data to be stored: number of edges of the eVIG
extern double *file_n_iter_ls;								// data to be stored: number of iterations of Local Search
extern double *file_n_best_fitness_improv;					// data to be stored: number of improvements of the best solution
extern double *file_n_escapes;								// data to be stored: number of escapes from local optimum
extern double *file_HD_LO;									// data to be stored: Hamming distance between two local optima
extern double *file_HD_pert;								// data to be stored: Hamming distance between solution before and after perturbation
extern double *file_fl_pert;								// data to be stored: number of subfunctions changed by perturbation
// Matrices
extern int **File_best_ind;									// data to be stored: best individual
	
// Function declaration
//  aux_functions.cpp
int *aloc_vectori(int lines);
double *aloc_vectord(int lines);
int **aloc_matrixi(int lines , int collums);
double **aloc_matrixd(int lines , int collums);
void desaloc_matrixi(int **Matrix , int lines);
void desaloc_matrixd(double **Matrix , int lines);
int random_int(int L_range, int H_range);
double random_dou(void);
void rand_perm_size(int *inp, int *out, int size_inp, int size_out);
// file_man.cpp
void file_output(int N, int K, int model_NK, int total_runs, int rball_size);
void read_problem(int N, int K, int instance, list<int> *PHI, double **Fl);
// statistics.cpp
void statistics(individual *sol, int n_run, long int n_iter_ls);
// selection.cpp
void AcceptanceCriterion(individual *x, individual *x_new);
// perturbation.cpp
void randomKickbP(int *x, int *x_new, int rball_size);

