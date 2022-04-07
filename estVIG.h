/******************************************************************************\
*				Class: estimated Variable Interaction Graph				 	   *
\******************************************************************************/
#include <cmath>
#include <cstdlib>

class estVIG {
	private:	
		void randPerm(int *out, int size_inp, int size_out);					// Get size_out elements of random Permutation of integers from 0 to size_inp-1	
		int randomInt(int L_range, int H_range);			 					// Random integer between L_range and H_range
	public:
	    int N;																	// Number of decision variables
	    list<int> *eVIG; 														// Vector of lists eVIG - estimated Variable Interaction Graph (VIG) (eVIG[i]: indicates the variables that interact with the i-th variable )	
		estVIG(int N);
	    ~estVIG(void);
	    void print(void);														// Print eVIG
	    void save(int N, int K, int model_NK, int n_run, int perturbation_type, int rball_size);						// Save eVIG
	    int addEdge(int a, int b);												// Add edge (a,b) to the eVIG
		int eVIGbP(int *x, int *x_new);											// Perturbation based on flipping bits of neighboring variables based on the eVIG
};


/******************************************************************************\
*								Constructor									   *
\******************************************************************************/
estVIG::estVIG(int N){	
	
	// Parameters for eVIG
	this->N = N;
						
	// Vector of lists eVIG - estimated Variable Interaction Graph (VIG) (eVIG[i]: indicates the variables that interact with the i-th variable )	
	eVIG = new list<int>[N]; 				
	
}


/******************************************************************************\
*								 Destructor									   *
\******************************************************************************/
estVIG::~estVIG(void){
	
	delete [] eVIG; 
	          
}  


/******************************************************************************\
*						 Random Integer between L_range and H_range			   *
\******************************************************************************/
int estVIG::randomInt(int L_range, int H_range)
{
	return(  (int) ( (rand()/(RAND_MAX+1.0))*(H_range-L_range+1)+L_range ) );  	// random integer beteween [L_range and H_range]
}


/******************************************************************************\
* Get size_out elements of random Permutation of integers from 0 to size_inp-1 *
\******************************************************************************/
void estVIG::randPerm(int *out, int size_inp, int size_out)
{
	int j, aux, *auxv;
	
	if (size_inp<size_out){
		cout<<"Error in function randPerm!"<<endl;
		exit(1);
	}
	
	auxv=new int[size_inp];
	
	for(int i=0;i<size_inp;i++)
		auxv[i]=i;	

	for(int i=0;i<size_out;i++) {
		j= randomInt(i,size_inp-1);  		
		aux=auxv[i];
		auxv[i]=auxv[j];
		auxv[j]=aux;
		out[i]=auxv[i];
	}
	
	delete [] auxv;
}


/*********************************************************************************\
* Perturbation based on flipping bits of neighboring variables based on the eVIG  *
\*********************************************************************************/
int estVIG::eVIGbP(int *x, int *x_new){
	int i, j, r, n_el=0, n_el_changed, max_n_el_changed, *order_v;

	for (i=0;i<N;i++)
		x_new[i]=x[i];
	
	// Randomly selecting first variable
	r=randomInt(0,N-1); 
	// transform x_new[r]
	if (x_new[r]==0)
		x_new[r]=1;
	else
		x_new[r]=0;
	n_el++;	
	if (eVIG[r].empty()){
		return n_el;
	}	
		
	// limiting the maximum number of variables that will be changed
	n_el_changed=eVIG[r].size();	
	max_n_el_changed=0.5*N-1;				// limiting the total number in 50% of N
	if (n_el_changed>max_n_el_changed)
		n_el_changed=max_n_el_changed;
	
	// Defining the order of the variables
	order_v=new int[eVIG[r].size()];	
	if (n_el_changed>0){		
		i=0;
		for(list<int>::iterator jj = eVIG[r].begin(); jj !=eVIG[r].end(); jj++){
			order_v[i]=*jj;
			i++;
		}	
		if (n_el_changed==max_n_el_changed){
			randPerm(order_v, eVIG[r].size(), n_el_changed);		// random permutation
		}
	}
	
	// Flipping neighboring variables
	if (n_el_changed>0){
		for (i=0;i<n_el_changed;i++){		
			j=order_v[i];
			// transform x_new[j]
			if (x_new[j]==0)
				x_new[j]=1;
			else
				x_new[j]=0;
			n_el++;
		}			
	}
	
	delete [] order_v;
	
	return n_el;			
}


/******************************************************************************\
*								Print Mk information						   *
\******************************************************************************/
void estVIG::print(void){

	cout<<"eVIG: variables "<<endl;
	for (int i=0;i<N;i++){
		cout<<" x_"<<i<<": ";
		for(list<int>::iterator j = eVIG[i].begin(); j != eVIG[i].end(); j++) 
			cout<<"  x_"<<*j;
		cout<<endl;
	}
		
}


/******************************************************************************\
*								Save NNet information														   *
\******************************************************************************/
void estVIG::save(int N, int K, int model_NK, int n_run, int perturbation_type, int rball_size){
	FILE *eVIG_file;
	char *name_p;
	char name[CHAR_LEN];

	//cout<<"Saving eVIG"<<endl;
    name_p = name;
	sprintf(name,"eVIG_N%d_K%d_m%d_p%d_s%d_r%d.dat",N,K,model_NK,perturbation_type,rball_size,n_run);
	if ((eVIG_file = fopen(name_p,"w"))==NULL) {
		puts("The file eVIG to be saved cannot be open \n");
		exit(1);
	}
		
	for (int i=0;i<N;i++){
		fprintf(eVIG_file," x_%d",i);
		for(list<int>::iterator j = eVIG[i].begin(); j != eVIG[i].end(); j++) 
			fprintf(eVIG_file," x_%d",*j);
		fprintf(eVIG_file,"\n");
	}

	fclose(eVIG_file);
		
}


/******************************************************************************\
* 				Add edge (a,b) to the eVIG									   *
\******************************************************************************/
int estVIG::addEdge(int a, int b){
	int aux_flag=1;
	list<int>::iterator ii;
	
	// Check if link (a,b) is already in the eVIG	
	if ( !eVIG[a].empty() && !eVIG[b].empty() ){
		ii=eVIG[a].begin();
		while (ii != eVIG[a].end() && *ii != b)					
			ii++;
		if (ii != eVIG[a].end())	
			aux_flag=0;
	}
	if (aux_flag==1){
		// Add edge (a,b)
		eVIG[a].push_back(b); 
		eVIG[b].push_back(a); 	
	}
	
	return aux_flag;
	
}
