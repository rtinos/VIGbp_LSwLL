/******************************************************************************\
*				  				 Files Manipulation							 *
\******************************************************************************/

#include "defs.h"
#include<cstring>
#include<fstream>


/*******************************************************************\
*	 Read the problem instance  								    *
* Obs.: the order (big endian) for the variables in the 		    *
*		file is different from the order (little endian) used here  *
\*******************************************************************/
void read_problem(int N, int K, int instance, list<int> *PHI, double **Fl){
	  int k=0, n_sub, *v_aux;	  
	  char line[CHAR_LEN], *keywords, Delimiters[] = " :=\n\t\r\f\v", name[CHAR_LEN], *name_p;
		
	  v_aux=aloc_vectori(K+1);
	  n_sub=(int) pow(2.0,K+1);
	  name_p = name;
	  sprintf(name,"prob/N%d/K%d/solver_fitcont_%d.dat",N,K,instance);  // files are in directory prob
	
	  ifstream fin(name_p);
	  while((fin.getline(line, CHAR_LEN-1))){
			if(!(keywords = strtok(line, Delimiters)))
	  			continue;
			if(!strcmp(keywords, "m")){
				for(int i=0; i<K+1; i++)
					v_aux[i] = atoi(strtok(NULL, Delimiters));
				for(int i=0; i<K+1; i++)
					PHI[k].push_back(v_aux[K-i]);				
	  			for(int i=0; i<n_sub; i++)
					fin>>Fl[k][i];
				k++;
			}
			
	  }
	  fin.close();
	  
	  delete [] v_aux;
}


/******************************************************************************\
* 					Save data : end of the simulation						   *
\******************************************************************************/
void file_output(int N, int K, int model_NK, int total_runs, int rball_size)
{
	char *name_p;
	char name[CHAR_LEN];
	FILE *Bestfit_file, *Bestind_file, *Time_file, *Gen_file, *Nit_file, *Imp_file;

    name_p = name;
		
  	// Save extra information
  	if (save_extra_info_flag==1){
  		
  		// Mean Best fitness in each generation 
  		FILE *Bfg_file;
		sprintf(name,"bfg_N%d_K%d_m%d_p%d_s%d.dat",N,K,model_NK,perturbation_type,rball_size);
		if ((Bfg_file = fopen(name_p,"w"))==NULL) {
			puts("The file bfg to be saved cannot be open \n");
			exit(1);
		}
		for (int i=0;i<max_gen;i++) {
			fprintf(Bfg_file,"%.14f ",file_best_fitness_gen[i]/total_runs);
		}
		fclose(Bfg_file);			
		
		// Mean Hamming distance between consecutive local optima	
		double gen_aux;
		if (gen<max_gen)
			gen_aux=gen-1.0;
		else
			gen_aux=max_gen-1.0;
		FILE *Hdlo_file;
		sprintf(name,"hdlo_N%d_K%d_m%d_p%d_s%d.dat",N,K,model_NK,perturbation_type,rball_size);
		if ((Hdlo_file = fopen(name_p,"w"))==NULL) {
			puts("The file hdlo to be saved cannot be open \n");
			exit(1);
		}
		for (int i=0;i<total_runs;i++) {
				if (file_n_escapes[i]>0)
					fprintf(Hdlo_file,"%.14f ",file_HD_LO[i]/file_n_escapes[i]);
				else
					fprintf(Hdlo_file,"%.14f ",0.0);
		}
		fclose(Hdlo_file);
								
		// Mean number of escapes from local optimum
		FILE *Nesc_file;
		sprintf(name,"nesc_N%d_K%d_m%d_p%d_s%d.dat",N,K,model_NK,perturbation_type,rball_size);
		if ((Nesc_file = fopen(name_p,"w"))==NULL) {
			puts("The file nesc to be saved cannot be open \n");
			exit(1);
		}
		for (int i=0;i<total_runs;i++) {
				fprintf(Nesc_file,"%.14f ",file_n_escapes[i]/gen_aux);
		}
		fclose(Nesc_file);
		
		// Mean Hamming distance between solution before and after perturbation
		FILE *Hdpert_file;
		sprintf(name,"hdpert_N%d_K%d_m%d_p%d_s%d.dat",N,K,model_NK,perturbation_type,rball_size);
		if ((Hdpert_file = fopen(name_p,"w"))==NULL) {
			puts("The file hdpert to be saved cannot be open \n");
			exit(1);
		}
		for (int i=0;i<total_runs;i++) {
				fprintf(Hdpert_file,"%.14f ",file_HD_pert[i]/gen_aux);
		}
		fclose(Hdpert_file);
		
		// Mean number of subfunctions changed by perturbation
		FILE *Flpert_file;
		sprintf(name,"flpert_N%d_K%d_m%d_p%d_s%d.dat",N,K,model_NK,perturbation_type,rball_size);
		if ((Flpert_file = fopen(name_p,"w"))==NULL) {
			puts("The file flpert to be saved cannot be open \n");
			exit(1);
		}
		for (int i=0;i<total_runs;i++) {
				fprintf(Flpert_file,"%.14f ",file_fl_pert[i]/gen_aux);
		}
		fclose(Flpert_file);
	
	
	}
	if (perturbation_type==2){
		// Mean number of edges of the eVIG
		FILE *Nedges_file;
		sprintf(name,"nedges_N%d_K%d_m%d_p%d_s%d.dat",N,K,model_NK,perturbation_type,rball_size);
		if ((Nedges_file = fopen(name_p,"w"))==NULL) {
			puts("The file nedges to be saved cannot be open \n");
			exit(1);
		}
		for (int i=0;i<max_gen;i++) {					
			fprintf(Nedges_file,"%.3f ",file_n_edges_eVIG_gen[i]/total_runs);
		}		
		fclose(Nedges_file);
	}
	


    // Best fitness 
	sprintf(name,"bfi_N%d_K%d_m%d_p%d_s%d.dat",N,K,model_NK,perturbation_type,rball_size);
	if ((Bestfit_file = fopen(name_p,"w"))==NULL) {
		puts("The file bfi to be saved cannot be open \n");
		exit(1);
	}
	for (int i=0;i<total_runs;i++) {
		fprintf(Bestfit_file,"%.14f ",file_best_fitness[i]);
	}
	fclose(Bestfit_file);
		
	 // Best individuals
	sprintf(name,"bind_N%d_K%d_m%d_p%d_s%d.dat",N,K,model_NK,perturbation_type,rball_size);
	if ((Bestind_file = fopen(name_p,"w"))==NULL) {
		puts("The file bind to be saved cannot be open \n");
		exit(1);
	}
	for (int i=0;i<total_runs;i++) {
		for (int gene=0;gene<chrom_size;gene++)
			fprintf(Bestind_file,"%d ",File_best_ind[i][gene]);
		fprintf(Bestind_file,"\n");
	}
	fclose(Bestind_file);
	
  	// Time for each run
	sprintf(name,"time_N%d_K%d_m%d_p%d_s%d.dat",N,K,model_NK,perturbation_type,rball_size);
	if ((Time_file = fopen(name_p,"w"))==NULL) {
		puts("The file time to be saved cannot be open \n");
		exit(1);
	}
	for (int i=0;i<total_runs;i++) {
		fprintf(Time_file,"%.2f ",time_run[i]);
	}
	fclose(Time_file);
	
	// Number of generations for each run
	sprintf(name,"gen_N%d_K%d_m%d_p%d_s%d.dat",N,K,model_NK,perturbation_type,rball_size);
	if ((Gen_file = fopen(name_p,"w"))==NULL) {
		puts("The file gen to be saved cannot be open \n");
		exit(1);
	}
	for (int i=0;i<total_runs;i++) {
		fprintf(Gen_file,"%d ",file_gen[i]);
	}
	fclose(Gen_file);		
	
	// Number of iterations of Local Search
	sprintf(name,"nit_N%d_K%d_m%d_p%d_s%d.dat",N,K,model_NK,perturbation_type,rball_size);
	if ((Nit_file = fopen(name_p,"w"))==NULL) {
		puts("The file nit to be saved cannot be open \n");
		exit(1);
	}
	for (int i=0;i<total_runs;i++) {
		fprintf(Nit_file,"%.2f ",file_n_iter_ls[i]);
	}
	fclose(Nit_file);	

	// Number of improvements of the best solution
	sprintf(name,"imp_N%d_K%d_m%d_p%d_s%d.dat",N,K,model_NK,perturbation_type,rball_size);
	if ((Imp_file = fopen(name_p,"w"))==NULL) {
		puts("The file imp to be saved cannot be open \n");
		exit(1);
	}
	for (int i=0;i<total_runs;i++) {
		fprintf(Imp_file,"%.5f ",file_n_best_fitness_improv[i]);
	}
	fclose(Imp_file);	
	
}


