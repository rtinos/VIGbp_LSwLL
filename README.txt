*** Iterated Local Search (ILS) with Perturbation based on the Variable Interaction Graph (VIGbp) ***

Description: This is the source code for running ILS with VIGbp in the black box experiments with NK landscapes described in the paper:

Reference:  Tinos, R.; Przewozniczek, M. D.; Whitley, D. (2022), "Iterated Local Search with Perturbation based on Variables Interaction for Pseudo-Boolean Optimization'', Proc. of the Genetic and Evolutionary Computation Conference (GECCO '22), https://doi.org/10.1145/3512290.3528716.      

It also includes the code for Local Search with Linkage Learning (LSwLL).

Contact: Renato Tinos <rtinos@ffclrp.usp.br>


Running the code: ./ILS_NKland  <N> <K> <model_NK> <perturbation_type> <rball_size> (only for perturbation types 1, 4, and 5)

<N>: size of the NK landscapes instance. N is a positive integer.

<K>: controls the epistasis degree of the NK landscapes instance. K is a non-negative integer smaller than or equal to N.

<model_NK>: NK landscapes model. 0-random NK; 1-regular NK (in this case, a table with values of subfuncions f_l for instances solved by dynamic programming is loaded); 2-random NKq with q=2^(K+1)/2-1; ; 3-random NKq with q=2^(K+1)-1.

<perturbation_type>: type for perturbation in ILS. 1. SRP: standard random perturbation with fixed strength (rball_size); 2. estimated VIG based perturbation (black box VIGbp); 3. real VIG based perturbation (gray box VIGbp); 4. RPR: random perturbation with random strength.; 5. RPN: random perturbation with strength dependent on dimension N.

<rball_size> (only for perturbation types 1, 4, and 5): perturbation strength (rball_size)


Example for running the code for: random NK landscapes with N=100, K=2, black box VIGbp

make

./ILS_NKland 100 2 0 2


Observation 1: VIGb for Gray Box  is implemented in Class Mk, given in Mk.h. MK.h implements Mk landscapes (that is a generalization of NK landscapes). Some methods in Mk.h:

- double Mk::compFitness (int *x): Fitness computation for Mk landscapes. The version here does not use information of the VIG to efficiently compute the fitness, i.e., the fitness is computed as in Black Box problems.

- void Mk::createVIG(void): Create the (real) Variable Interaction Graph based on the informations of Mk landscapes.

- void Mk::VIGbP(int *x, int *x_new): gray box VIGbp, i.e., implements perturbation based on the real VIG. 


Observation 2: VIGb for Black Box is implemented in Class estVIG, given in estVIG.h. Some methods in estVIG.h:

- int estVIG::addEdge(int a, int b): Add edge (a,b) to the estimated VIG (eVIG).

- void estVIG::eVIGbP(int *x, int *x_new): black box VIGbp, i.e., implements perturbation based on the estimated VIG. 

- void estVIG::save(int N, int K, int model_NK, int n_run, int perturbation_type, int rball_size): save the estimated VIG.


Observation 3: Local Search (LS) and Local Search with Linkage Learning (LSwLL) are implemented in ILS.cpp:

- double LS(int *x, double f, long int *n_iter, Mk *Mk_instance): Standard (First Improvement) Local Search.

- double LSwLL(int *x, double f, long int *n_iter, Mk *Mk_instance, estVIG *eVIG_instance, int *n_edges_eVIG): (First Improvement) Local Search with Linkage Learning.


Observation 4: The parameters of the experiment (number of instances, number of runs for each instance) are defined in global.cpp.


Observation 5: ILS_NKland generates a file for each measure described in the paper. For example:
 
- bfi_N%d_K%d_m%d_p%d_s%d.dat (example for ./ILS_NKland 100 2 0 2: , bfi_N100_K2_m0_p2_s0.dat): save the best fitness found in each run

