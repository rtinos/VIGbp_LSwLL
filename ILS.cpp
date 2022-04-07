/*******************************************************************************************\
*  	Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 		  	   				*
*																							*
* 	Iterated Local Search with Perturbation based on the Variable Interaction Graph (VIG)	*
* 	Test problem: NK Landscapes						 										*
*						 																	*
* 	Copyright (C) 2022  Renato Tinos <rtinos@ffclrp.usp.br>									*
* 						 																	*
* 	Reference:  Tinos, R.; Przewozniczek, M. D.; Whitley, D. (2022), ``Iterated Local 		*
				Search with Perturbation based on Variables Interaction for Pseudo-Boolean 	*
				Optimization'', Proc. of the Genetic and Evolutionary Computation 		    *
				Conference (GECCO '22), https://doi.org/10.1145/3512290.3528716.           	*                    
* 						 																	*
* 	eVIGbpILS is free software: you can redistribute it and/or modify it 					*
* 		under the terms of the GNU General Public License as published by the				*
* 		Free Software Foundation, either version 3 of the License, or						*
* 		(at your option) any later version.						 							*
* 						 																	*
* 	eVIGbpILS is distributed in the hope that it will be useful, but						*
* 		WITHOUT ANY WARRANTY; without even the implied warranty of							*
* 		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.								*
* 		See the GNU General Public License for more details.								*
* 																							*
* 	You should have received a copy of the GNU General Public License along					*
* 		with this program.  If not, see <http://www.gnu.org/licenses/>.						*		
\*******************************************************************************************/

#include <time.h>
#include "defs.h"
#include "Mk.h"							// Mk landscapes class
#include "estVIG.h"						// estimated VIG class


/******************************************************************************\
*				  	Print population					 			 			*
\******************************************************************************/
void print_data(individual *sol, int n_run){

	cout <<"Generation:"<< gen << ", run: "<<n_run<<endl;
	cout <<"Fitness of the best solution:"<< file_best_fitness[n_run] << endl;
	cout <<"("<< sol->fitness<<") " ;
	/*for (int gene=0;gene<chrom_size ;gene++) 
		cout << sol->solutionVector[gene]<<" ";
	cout << endl;*/
	//system("pause");
	
}


/******************************************************************************\
*		Record Information - Statistics about optimization					   *
\******************************************************************************/
void computeStatOpt(Mk *Mk_instance, individual *sol, individual *sol_afterPert, individual *sol_afterLS, int n_run)
{
	int *v_move;
	double n_changed_bits;
	
	// Statistics for changing from sol to sol_afterLS
	n_changed_bits=0.0;
	for (int gene=0;gene<chrom_size;gene++){
		if (sol->solutionVector[gene]!=sol_afterLS->solutionVector[gene])
				n_changed_bits+=1.0;
	}
	if (n_changed_bits>0){
		file_n_escapes[n_run]=file_n_escapes[n_run]+1.0;
		file_HD_LO[n_run]=file_HD_LO[n_run]+n_changed_bits;
	}
	
	// Statistics for changing from sol to sol_afterPert
	v_move=aloc_vectori(chrom_size);
	n_changed_bits=0.0;
	for (int gene=0;gene<chrom_size;gene++){		
		if (sol->solutionVector[gene]!=sol_afterPert->solutionVector[gene]){
				n_changed_bits+=1.0;
				v_move[gene]=1;
		}
		else
			v_move[gene]=0;		
	}
	if (n_changed_bits>0){
		file_HD_pert[n_run]=file_HD_pert[n_run]+n_changed_bits;
		file_fl_pert[n_run]=file_fl_pert[n_run]+ ((double) Mk_instance->compFlChange(v_move));
	}
	delete [] v_move;
	
}


/******************************************************************************\
*								Fitness (real) Computation    				   *
\******************************************************************************/
double compFitness( Mk *Mk_instance, variableType *x ){
		double Fitness;

		Fitness = Mk_instance->compFitness(x);		// fitness function for Mk landscapes

		return Fitness;
}


/******************************************************************************\
*					Local Search: First Improvement							   *
\******************************************************************************/
double LS(int *x, double f, long int *n_iter, Mk *Mk_instance){
	int i, j=0, ni=0, nr=0, *R;
	double df, fit_aux;

	// Randomly defining the order of the search
	R=aloc_vectori(chrom_size);							// random sequence for visiting variables
	rand_perm_size(vsort_aux, R, chrom_size, chrom_size);
	
	// Local Search	
	*n_iter=0;
	while(ni<chrom_size){
		i=R[j];
		// Flipping the i-th bit
		if (x[i]==0)
			x[i]=1;
		else
			x[i]=0;
		
		// Computing the difference in the fitness
		fit_aux=compFitness(Mk_instance,x);			// fitness computation
		df=fit_aux-f;
		ni++;
		
		// Cheking improvement
		if (df<EPS1){	
			//  No Improvement: Reversing the flip
			if (x[i]==0)
				x[i]=1;
			else
				x[i]=0;
		}
		else{
			// Improvement
			f=fit_aux;
			ni=0;
			nr++;
			if (nr>5*chrom_size)
				ni=chrom_size;					// a limit for the number of repetitions (improvements) is here defined					
		}
		j++;
		if (j==chrom_size)
			j=0;
		*n_iter=*n_iter+1;
										 	
	}
	
	delete [] R;
	
	return f;

}


/******************************************************************************\
*		Local Search with Linkage Learning (LSwLL): First Improvement		   *
\******************************************************************************/
double LSwLL(int *x, double f, long int *n_iter, Mk *Mk_instance, estVIG *eVIG_instance, int *n_edges_eVIG){
	int ni=0, nr=0, *R, g, h, *Q, delta, k_R, j_Q, r=0;
	double df, fit_aux, *F;

	// Randomly defining the order of the search
	R=aloc_vectori(chrom_size);								// R defines a random sequence for variables that were not visited yet
	rand_perm_size(vsort_aux, R, chrom_size, chrom_size);	// generate random sequence for variables that were not visited yet

	// Allocation of lists Q and F
	Q=aloc_vectori(chrom_size);								// list Q record variables already visited with no improvement
	F=aloc_vectord(chrom_size);								// list F records the differences in the fitness for variable already visited with no improvement

	// Initiating markers 
	k_R=0;													// k_R is the index of the current element of R; start at first position
	delta=0;												// delta is the number of variables already visited
	j_Q=0;													// j_Q is the index of the current element of Q; start at first position
	
	// Local Search	
	*n_iter=0;												// n_iter records number of iterations of LS
	while(ni<chrom_size){									// stop criterion
		// Defining g
		if (j_Q<delta){
			g=Q[j_Q];										// g controls the variable that will be changed; here, g receives a variable already visited
		}
		else{
			g=R[k_R];										// here, g receives a variable not visited yet
			k_R++;
			if (k_R==chrom_size)
				k_R=0;
		}
		
		// Flipping the g-th bit
		if (x[g]==0)
			x[g]=1;
		else
			x[g]=0;					
		// Computing the difference in the fitness
		fit_aux=compFitness(Mk_instance,x);					// fitness computation
		df=fit_aux-f;										// difference in the fitness caused by flipping g	
		ni++;												// ni is the number of flips before improvement 
		
		// Adding edges to eVIG when j_Q<delta
		if (j_Q<delta){										// variable is already visited 				
			if ( fabs(df-F[j_Q])>EPS2 ){					
				// df changed because variables g and h interact
				if (eVIG_instance->addEdge(g,h)){			// adding edge (g,h) to eVIG (if the edges does not exist)
					*n_edges_eVIG=*n_edges_eVIG+1;			// update number of edges in eVIG
				}			
			}
		}
					
		// Cheking improvement
		if (df<EPS1){	
			// No Improvement
			// Reversing the flip 
			if (x[g]==0)
				x[g]=1;
			else
				x[g]=0;	
								
			// update Q
			if (j_Q>=delta){
				// variable is not visited yet
				Q[j_Q]=g;								// add new visited variable to list Q
				F[j_Q]=df;							// records its df				 	
			}			
			j_Q++;										// update number of visited variables										
		}
		else{
			// Improvement				
			f=fit_aux;		
			h=g;										// record variable with improvement			
			if (j_Q<delta)
				delta=0;								// if the improvement occurred in a visited variable, delta is set to zero, i.e., Q and F are emptied					
			else	
				delta=j_Q;								// if the improvement occurred in a non-visited variable, delta is set to j_Q
				
			// revisit sequence only one time (r is the visiting flag)
			if (r==1){
				delta=0;
				r=0;
			}
			else{
				r=1;
			}
			
			j_Q=0;								
			ni=0;								// ni (the number of flips before improvement) is set to zero
			nr++;								// nr is the number of improvements
			if (nr>5*chrom_size)
				ni=chrom_size;					// a limit for the number of repetitions (improvements) is here defined					
		}		
				
		*n_iter=*n_iter+1;
										 	
	}
	
	delete [] R;
	delete [] Q;
	delete [] F;
	
	return f;

}


/******************************************************************************\
*				  	Perturbation				   		 				 	   *
\******************************************************************************/
void Perturbation(individual *sol, individual *sol_new, Mk *Mk_instance, int rball_size, long int n_iter_ls, estVIG *eVIG_instance){

	// Generating new solution
	if ( perturbation_type==1 ){
		// Perturbation based on randomly kicking rball_size bits (SRP)	
		randomKickbP(sol->solutionVector, sol_new->solutionVector, rball_size);
	}
	else if ( perturbation_type==2 ){
		// Perturbation based on flipping bits of neighboring variables based on the estimated VIG (black box VIGbp)
		if (eVIG_instance->eVIGbP(sol->solutionVector, sol_new->solutionVector)==1){
			// Perturbation based on kicking rball_size bits based on the VIG
			randomKickbP(sol->solutionVector, sol_new->solutionVector, random_int(1, 0.1*chrom_size));  // if randomly selected bit in eVIGbP does not have neigbours in eVIG, then apply random perturbation 
		}
	}
	else if ( perturbation_type==3 ){
		// Perturbation based on flipping bits of neighboring variables based on the real VIG (gray box VIGbp)
		Mk_instance->VIGbP(sol->solutionVector, sol_new->solutionVector);	
	}
	else if ( perturbation_type==4 ){
		// Perturbation based on randomly kicking 1<=random x<=rball_size bits (RPR) 
		randomKickbP(sol->solutionVector, sol_new->solutionVector, random_int(1, rball_size));	
	}
	else if ( perturbation_type==5 ){
		// Perturbation based on randomly kicking rball_size*chrom_size/100.0 bits	(RPN)
		int rball_size_aux;
		rball_size_aux=rball_size*chrom_size/100.0;
		if (rball_size_aux>chrom_size)
			rball_size_aux=chrom_size;
		else if (rball_size_aux<1)
			rball_size_aux=1;
		randomKickbP(sol->solutionVector, sol_new->solutionVector, rball_size_aux);				
	}
	
	// Evavaluating the fitness of the new solution
	sol_new->fitness=compFitness(Mk_instance,sol_new->solutionVector);
		
			
}


/******************************************************************************\
*				  	Generate Initial Solution   		 				 	   *
\******************************************************************************/
void GenerateInitialSolution(individual *sol, Mk *Mk_instance){
				
	// Random Initialization
	for (int gene=0;gene<chrom_size;gene++) 
     	sol->solutionVector[gene] = random_int(0,1);
    sol->fitness = compFitness(Mk_instance, sol->solutionVector);		
	/*cout <<"("<< sol->fitness<<") " ;
	for (int gene=0;gene<chrom_size ;gene++) 
		cout << sol->solutionVector[gene]<<" ";
	cout << endl;
	system("pause");  	*/
}


/******************************************************************************\
*				  	Run of the ILS			 								   *
* Reference: Lourenco et al (2019). "Iterated Local Search: Framework and	   *
*				Applications", Handbook of Metaheuristics. Springer, 129-168.  *
\******************************************************************************/
void ils(Mk *Mk_instance, int n_run, int rball_size, int K, int model_neig ){
	int n_edges_eVIG=0;
	long int n_iter_ls;
	double time_aux, max_time; 
	individual s_star, s_prime, s_star_prime;
	clock_t time_start;
				
	// Initialization
	estVIG *eVIG_instance = new estVIG(chrom_size);				// from class estVIG (estVIG.h)
	if (perturbation_type==2){	
		Mk_instance->createVIG();								// here, VIG is created only for comparison purpose
		// Number of edges in VIG
		total_edges_VIG=Mk_instance->n_edges_VIG();	
	}
	max_time=((double) (Mk_instance->M*Mk_instance->k)/40.0);	
	time_start=clock();		
	s_star.solutionVector = aloc_vectori(chrom_size);
	s_prime.solutionVector = aloc_vectori(chrom_size);	
	s_star_prime.solutionVector = aloc_vectori(chrom_size);
	if (save_extra_info_flag==1){
		file_n_escapes[n_run]=0.0;
		file_HD_LO[n_run]=0.0;
		file_HD_pert[n_run]=0.0;
		file_fl_pert[n_run]=0.0;
	}
	file_n_iter_ls[n_run]=0.0;
	file_n_best_fitness_improv[n_run]=0.0;	
	gen=0;
	if (perturbation_type==3)	
		Mk_instance->createVIG();
	//Mk_instance->print();											// print graph and costs of the elements for the NK Landscape

	// Generate initial solutions				
	GenerateInitialSolution(&s_star, Mk_instance);															// generate initial solution
	//cout<<"Initialization: Generation: "<<gen<<", fitness before local search:"<<s_star.fitness<<endl;
	if (perturbation_type==2){
		s_star.fitness = LSwLL(s_star.solutionVector, s_star.fitness, &n_iter_ls, Mk_instance,eVIG_instance, &n_edges_eVIG); // local search with Linkage Learning: first improvement 	
		if (n_run==0)
			file_n_edges_eVIG_gen[gen]=( (double) n_edges_eVIG)/total_edges_VIG;
		else 
			file_n_edges_eVIG_gen[gen]+=( (double) n_edges_eVIG)/total_edges_VIG;
	}										
	else{
		s_star.fitness = LS(s_star.solutionVector, s_star.fitness, &n_iter_ls, Mk_instance); 						// local search: first improvement 				
	}
	
	//cout<<"...after "<<n_iter_ls<<" iterations of local search: "<<s_star.fitness<<endl;
	statistics(&s_star,n_run,n_iter_ls);
	//print_data(&s_star,n_run);
	
	// Generate new solutions	
	do {
		gen++; 																								// generation index
		
		Perturbation(&s_star, &s_prime, Mk_instance, rball_size, n_iter_ls, eVIG_instance);
		for (int gene=0;gene<chrom_size;gene++) 
     		s_star_prime.solutionVector[gene] = s_prime.solutionVector[gene];
		//cout<<"Perturbation: Generation: "<<gen<<", fitness before local search:"<<s_prime.fitness<<endl;
		if (perturbation_type==2){
			s_star_prime.fitness = LSwLL(s_star_prime.solutionVector, s_prime.fitness, &n_iter_ls, Mk_instance, eVIG_instance, &n_edges_eVIG); 	 // local search with Linkage Learning: first improvement (for Mk landscapes)			
			if (n_run==0)
				file_n_edges_eVIG_gen[gen]=( (double) n_edges_eVIG)/total_edges_VIG;
			else 
				file_n_edges_eVIG_gen[gen]+=( (double) n_edges_eVIG)/total_edges_VIG;
		}
		else{
			s_star_prime.fitness = LS(s_star_prime.solutionVector, s_prime.fitness, &n_iter_ls, Mk_instance); 				// local search: first improvement (for Mk landscapes)
		}
	
		if (save_extra_info_flag==1 && gen<max_gen)
			computeStatOpt(Mk_instance, &s_star, &s_prime, &s_star_prime, n_run);			// statistics about optimization
	
		AcceptanceCriterion(&s_star, &s_star_prime);
		
		statistics(&s_star,n_run,n_iter_ls);	
		//print_data(&s_star,n_run);					
		time_aux = ( (double) ( clock() - time_start ) ) / ( (double) CLOCKS_PER_SEC);

	//}while ( time_aux < max_time );		// for experiments with fixed time
	}while ( gen < max_gen );    			// for experiments with fixed number of iterations
		
	// Data to be saved
	time_run[n_run]=time_aux;
	file_gen[n_run]=gen;
	file_n_iter_ls[n_run] = file_n_iter_ls[n_run]/(gen+1.0);
	if (gen<max_gen)
		file_n_best_fitness_improv[n_run] = file_n_best_fitness_improv[n_run]/(gen-1.0);
	else
		file_n_best_fitness_improv[n_run] = file_n_best_fitness_improv[n_run]/(max_gen-1.0);
	
	
	if (perturbation_type==2){
		//eVIG_instance->print();
		// Save eVIG
		//eVIG_instance->save(chrom_size, K, model_neig, n_run, perturbation_type, rball_size);
		// Compare VIG and eVIG
		if (Mk_instance->compareVIG(eVIG_instance->eVIG)){
			exit(1);
		}
		
	}
			
	delete [] s_prime.solutionVector;
	delete [] s_star.solutionVector;
	delete [] s_star_prime.solutionVector;
	delete eVIG_instance;

}


/******************************************************************************\
*				  	NK model  												   *
\******************************************************************************/
void NK(int N, int K, int model_neig, int instance, list<int> *PHI, double **Fl) 
{
	int i, n_comb_fl, *v_aux, q_NK_1; 
	
	n_comb_fl=(int) pow(2.0,K+1);   					// Fl is a matrix N x n_comb_fl

	if (model_neig==1){
		// adjacent NK model
		read_problem(N, K, instance, PHI, Fl);			// read Fl and VIG for adjacent model file
	}
	else{
		if (model_neig==2 || model_neig==3){
			// random NKq model
			if (model_neig==2)
				q_NK_1=n_comb_fl/2;
			else
				q_NK_1=n_comb_fl;
			q_NK_1=q_NK_1-1;			
			// matrix Fl
			for (int l=0;l<N;l++)
				for (int j=0;j<n_comb_fl;j++)
					Fl[l][j]=(double) random_int(0, q_NK_1)/q_NK_1;	// random fl 
		}
		else{
			// random NK model
			// matrix Fl
			for (int l=0;l<N;l++)
				for (int j=0;j<n_comb_fl;j++)
					Fl[l][j]=random_dou();					// random fl 
		}
		
	    // vector of lists PHI (PHI[l] indicates the variables that influences subfunction fl)
		v_aux=aloc_vectori(K+1);
		for (int l=0;l<N;l++){
			PHI[l].push_back(l); 						// index l is always in PHI[l]
			rand_perm_size(vsort_aux, v_aux, N, K+1);	
			i=0;	
			for (int j=0;j<K;j++){
				if (v_aux[i]==l)
					i++;				
				PHI[l].push_back(v_aux[i]); 			// random index in PHI[l]
				i++;
			}
		}
		delete [] v_aux;
	}

}


/******************************************************************************\
*				  	Main													   *
\******************************************************************************/
int main(int argc , char *argv[])
{
	int total_runs, n_run=0, rball_size;
	int N_NK, K_NK, model_NK;					// variables for the NK model
	double **Fl_NK;								// 	matrix with the values for the subfunctions fl in the NK model
	list<int> *PHI_NK;							// 	vector of lists PHI (PHI[l] indicates the variables that influences subfunction fl)
	
	// Arguments
	if( argc < 5) {
		cout<<"Insufficient number of arguments!"<<endl;
		cout<<"Call: ILS_NKland <N> <K> <model_NK> <perturbation_type> ( for perturbation_type!=2 and perturbation_type!=3: <rball_size>)"<<endl;
		exit(1);
	}
	else{
		N_NK=atoi(argv[1]);
		K_NK=atoi(argv[2]);	
		model_NK=atoi(argv[3]);
		perturbation_type=atoi(argv[4]);							
		if ( N_NK<1 || K_NK<0 || K_NK>N_NK-1 || model_NK<0 || model_NK>3 || perturbation_type<1 || perturbation_type>5 ) {
			cout<<"Incorrect arguments!"<<endl;
			cout<<"Call: ILS_NKland  <N> (N>0) <K> (0<=K<=N) <model_NK> (0-random NK; 1-regular NK; 2,3-random NKq) <perturbation_type> (1<=perturbation_type<=5) ( for perturbation_type !=2 and !=3: <rball_size> (>0 && <N))"<<endl;
			exit(1);
		}
		if ( perturbation_type!=2 && perturbation_type!=3){ 			
			if( argc < 6) {
				cout<<"Insufficient number of arguments!"<<endl;
				cout<<"Call: ILS_NKland <N> <K> <model_NK> <perturbation_type> ( for perturbation_type!=2 and !=3: <rball_size>)"<<endl;
				exit(1);
			}
			rball_size=atoi(argv[5]);
			if ( rball_size<2)
				rball_size=2;
			else if(perturbation_type!=5 && rball_size>N_NK) 
				rball_size=N_NK;
			else if(perturbation_type==5 && rball_size>100) 
				rball_size=100;	
		}
		else{
			rball_size=0;
		}
	}	
	
	// Parameters
	chrom_size=N_NK;													// size of the solutionVector
	total_runs=n_instances*n_runs_per_instance;							// number of runs
	max_gen=5000;														// maximum number of generations (when stop criterion is time, it is used for saving some statistics data)	
			
	// Allocation of vectors and matrices
	Fl_NK = aloc_matrixd ( N_NK , (int) pow(2.0,K_NK+1) );
	file_best_fitness_gen=aloc_vectord(max_gen);
	file_n_edges_eVIG_gen=aloc_vectord(max_gen);
	file_best_fitness=aloc_vectord(total_runs);
	file_gen=aloc_vectori(total_runs);
	file_n_iter_ls=aloc_vectord(total_runs);
	file_n_best_fitness_improv=aloc_vectord(total_runs);
	time_run=aloc_vectord(total_runs);			
	file_n_escapes=aloc_vectord(total_runs);
	file_HD_LO=aloc_vectord(total_runs);
	file_HD_pert=aloc_vectord(total_runs);
	file_fl_pert=aloc_vectord(total_runs);
	File_best_ind=aloc_matrixi(total_runs,chrom_size);
	vsort_aux=aloc_vectori(chrom_size);									// Auxiliar sorted vector of integers (used in different methods)
	for (int i=0;i<chrom_size;i++)
		vsort_aux[i]=i;
	if (save_extra_info_flag==1){
		for (int i=0;i<max_gen;i++)
			file_best_fitness_gen[i]=0.0;
	}
	cout << "\n ***** Iterated Local Search ****" << endl;
	cout << "NK landscapes with: N="<<N_NK <<", K="<<K_NK<< endl;
	cout << "Perturbation Method="<<perturbation_type<<", rball_size="<<rball_size<< endl;

	for (int instance=0;instance<n_instances;instance++) {	
		cout <<"Instance: "<< instance << endl;
		srand(instance);												// random seed  (for instance)	
			
		// Creating the NK Model
		PHI_NK = new list<int>[N_NK]; 
		NK(N_NK, K_NK, model_NK, instance, PHI_NK, Fl_NK);				// create matrices VIG_NK and Fl_NK
		// Creating the Mk Model	
		Mk *Mk_instance = new Mk(N_NK, N_NK, K_NK+1, PHI_NK, Fl_NK);	// from class Mk (Mk.h)
		
		for (int i_run=0;i_run<n_runs_per_instance;i_run++){
			cout <<"Run:"<< n_run <<", random seed: " << n_instances+n_run << endl;
			srand( n_instances+n_run );									// random seed  (for run)
			ils(Mk_instance, n_run, rball_size, K_NK, model_NK);	    // run ILS
			n_run++;
		}

		delete [] PHI_NK;
		delete Mk_instance;
	}
	
	file_output(N_NK,K_NK,model_NK,total_runs, rball_size);					// save data

	// Desallocation of vectors and matrices	
	desaloc_matrixd (Fl_NK,N_NK);
	desaloc_matrixi (File_best_ind,total_runs);
	delete [] time_run;
	delete [] file_n_iter_ls;
	delete [] file_n_best_fitness_improv;
	delete [] file_gen;
	delete [] file_best_fitness;
	delete [] file_best_fitness_gen;
	delete [] file_n_edges_eVIG_gen;
	delete [] vsort_aux;
	delete [] file_n_escapes;
	delete [] file_HD_LO;

	return 0;
}
