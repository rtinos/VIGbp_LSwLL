/******************************************************************************\
*								 Class: Mk Landscapes					 	   *
\******************************************************************************/
#include <cmath>
#include <cstdlib>

class Mk {
	private:	
		list<int> *PHI; 														// Vector of lists PHI (PHI[l]: indicates the variables that influences subfunction fl)	
		list<int> *PSI; 														// Vector of lists PSI (PSI[i]: indicates the indices of the subfunction fl influenced by the i-th variable)	
		list<int> *VIG; 														// Vector of lists VIG (VIG[i]: indicates the variables that interact with i-th variable in subfunctions fl)	
		double **Fl;															// Matrix with the contributions for each subfunction fl
		void randPerm(int *out, int size_inp, int size_out);					// Get size_v elements of random Permutation of integers from 0 to size_inp-1		
		int binvec2dec(int *v, int size_v);										// Transform a binary vector into a decimal (big endian notation)
		int randomInt(int L_range, int H_range);			 					// Random integer between L_range and H_range
		double dfElement(int *x, int var_index);								// Fitness difference caused by changing bit
	public:
	    int N;																	// Number of decision variables
	    int M;																	// Number of subfunctions
		int k;																	// Epistasis degree
		Mk(int N, int M, int k, list<int> *PHI, double **Fl);
	    ~Mk(void);
	    void createVIG(void);													// Create the Variable Interaction Graph 
	    void print(void);														// Print Mk information
		int compFlChange(int *v);												// Compute number of subfunctions changed by move v
		int compareVIG(list<int> *VIG_n);										// Compare VIG
		int n_edges_VIG( );														// Total number of edges in the VIG
	    double compFitness(int *x);												// Compute Fitness
		void VIGbP(int *x, int *x_new);											// Perturbation based on flipping bits of neighboring variables based on the VIG 
};


/******************************************************************************\
*								Constructor									   *
\******************************************************************************/
Mk::Mk(int N, int M, int k, list<int> *PHI, double **Fl){
	int var_index;
	
	// Parameters for Mk landscape
	this->N = N;
	this->k = k;
	this->M = M;
				
	// Vector of lists PHI (PHI[l]: indicates the variables that influences subfunction fl)		
	this->PHI = PHI;

	// Vector of lists PSI (PSI[i]: indicates the index of the subfunction fl influenced by the i-th variable)
	PSI = new list<int>[N]; 
	for (int l=0;l<M;l++){
		for(list<int>::iterator j = PHI[l].begin(); j != PHI[l].end(); j++){
			var_index=*j;
			PSI[var_index].push_back(l);
		}
	}	
		
	// Vector of lists VIG (VIG[i]: indicates the variables that interact with i-th variable in subfunctions fl)
	VIG = new list<int>[N]; 
				
	// Matrix with the contributions for each f_l
	this->Fl = Fl;
		
}


/******************************************************************************\
*								 Destructor									   *
\******************************************************************************/
Mk::~Mk(void){
	
	delete [] VIG; 
	delete [] PSI;  
	          
}  


/******************************************************************************\
* Create the Variable Interaction Graph 									   *
\******************************************************************************/
void Mk::createVIG(void){
	int i, i_next, aux_flag, var_index, var_index_next, *v_aux;
	list<int>::iterator ii;
	
	// Vector of lists VIG (VIG[i]: indicates the variables that interact with i-th variable in subfunctions fl)
	v_aux=new int[k];
	for (int l=0;l<M;l++){
		i=0;
		for(list<int>::iterator j = PHI[l].begin(); j != PHI[l].end(); j++){
			v_aux[i]=*j;
			i++;
		}
		for (i=0;i<k;i++){
			var_index=v_aux[i];
			i_next=i+1;
			while (i_next<k){
				var_index_next=v_aux[i_next];
				aux_flag=0;
				if ( !VIG[var_index].empty() && !VIG[var_index_next].empty() ){
					ii=VIG[var_index].begin();
					while (ii != VIG[var_index].end() && *ii != var_index_next)					
						ii++;
					if (ii != VIG[var_index].end())	
						aux_flag=1;
				}
				if (aux_flag==0){
					VIG[var_index].push_back(var_index_next); 
					VIG[var_index_next].push_back(var_index); 	
				}
				i_next++;
			}
		}
	}
	delete [] v_aux;
	
}


/******************************************************************************\
* Transform a binary vector into a decimal 									   *
\******************************************************************************/
int Mk::binvec2dec(int *v, int size_v){
	int y=0, base=1;
	
	for (int i=0;i<size_v;i++){	
		y += v[i]*base;
		base=base*2;	
	}
	
	return( y );  
}


/******************************************************************************\
*							Compute Fitness									   *
\******************************************************************************/
double Mk::compFitness (int *x){
	int i, *v_aux;
	double f=0.0;

	v_aux = new int[k];

	for (int l=0;l<M;l++){
		i=0;
		for(list<int>::iterator j = PHI[l].begin(); j != PHI[l].end(); j++){	
			v_aux[i]=x[*j];
			i++;
		}
		f = f + Fl[l][binvec2dec(v_aux,k)];

	}

	f = f/N;

	delete [] v_aux;
	return f;
}


/******************************************************************************\
*						 Random Integer between L_range and H_range			   *
\******************************************************************************/
int Mk::randomInt(int L_range, int H_range)
{
	return(  (int) ( (rand()/(RAND_MAX+1.0))*(H_range-L_range+1)+L_range ) );  	// random integer beteween [L_range and H_range]
}


/******************************************************************************\
* Get size_out elements of random Permutation of integers from 0 to size_inp-1  *
\******************************************************************************/
void Mk::randPerm(int *out, int size_inp, int size_out)
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



/******************************************************************************\
*		Fitness difference caused by changing bit 							   *
\******************************************************************************/
double Mk::dfElement(int *x, int var_index){
	int i, l, i_sel, *v_aux;
	double df=0.0;

	if (!PSI[var_index].empty()){
		v_aux = new int[k];
		
		for(list<int>::iterator li = PSI[var_index].begin(); li != PSI[var_index].end(); li++){
			l=*li;
			i=0;
			for(list<int>::iterator j = PHI[l].begin(); j != PHI[l].end(); j++){	
				if (*j==var_index)
					i_sel=i;				// index of var_index in comp_v
				v_aux[i]=x[*j];
				i++;
			}
			df-=Fl[l][binvec2dec(v_aux,k)];			// removing contribution
			// Flipping the bit
			if (v_aux[i_sel]==0)
				v_aux[i_sel]=1;
			else
				v_aux[i_sel]=0;
			df+=Fl[l][binvec2dec(v_aux,k)];			// adding contribution	
		}
		
		delete [] v_aux;
	}

	return df;
}




/*********************************************************************************\
* Perturbation based on flipping bits of neighboring variables based on the VIG   *
\*********************************************************************************/
void Mk::VIGbP(int *x, int *x_new){
	int i, j, r, n_el=0, n_el_changed, max_n_el_changed, *order_v;
	//double df;

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
	if (VIG[r].empty()){
		return;
	}	
	
	// limiting the maximum number of variables that will be changed
	n_el_changed=VIG[r].size();	
	max_n_el_changed=0.5*N-1;				// limiting the total number in 50% of N
	if (n_el_changed>max_n_el_changed)
		n_el_changed=max_n_el_changed;
	
	// Defining the order of the variables
	order_v=new int[VIG[r].size()];	
	if (n_el_changed>0){		
		i=0;
		for(list<int>::iterator jj = VIG[r].begin(); jj !=VIG[r].end(); jj++){
			order_v[i]=*jj;
			i++;
		}	
		if (n_el_changed==max_n_el_changed){
			randPerm(order_v, VIG[r].size(), n_el_changed);		// random permutation
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
		

	return;
			
}


/******************************************************************************\
*			Compute number of subfunctions changed by move v				   *
\******************************************************************************/
int Mk::compFlChange(int *v){
	int l, *fl_aux, nfl=0;

	fl_aux = new int[M];
	
	for (l=0;l<M;l++)
		fl_aux[l]=0;
		
	for (int var_index=0;var_index<N;var_index++){
		if (v[var_index]==1 && !PSI[var_index].empty()){	
			for(list<int>::iterator li = PSI[var_index].begin(); li != PSI[var_index].end(); li++){
				l=*li;
				fl_aux[l]=1;		// fl is changed
			}
		}
	}
	
	for (l=0;l<M;l++){
		nfl+=fl_aux[l];
		
	}

	
	delete [] fl_aux;
	return nfl;
}


/******************************************************************************\
*								Compare VIG									   *
\******************************************************************************/
int Mk::compareVIG(list<int> *VIG_n){
	int i, flag_dif=0;
	list<int>::iterator j, k;
	

	for (i=0;i<N && flag_dif==0;i++){
		j = VIG_n[i].begin();
		while (j != VIG_n[i].end()&& flag_dif==0) {
			k = VIG[i].begin();
			while (k != VIG[i].end() && *k != *j) {	
				k++;
			}	
			if (k == VIG[i].end())
				flag_dif=1;			
			j++;
		}
	}
	
	if (flag_dif==1){
		cout<<"There are at least one edge of VIG_n that is not present in VIG!"<<endl;
		cout<<"in VIG_n  (x_"<<i<<", x_"<<*j<<")"<<endl;
		cout<<"  x_"<<*j;
		cout<<"in VIG x_"<<i<<": ";
		for(k = VIG[i].begin(); k != VIG[i].end(); k++) 
			cout<<"  x_"<<*k;
		cout<<endl;
	}
	
	return flag_dif;   		// return 1 if VIG_n has an edge not present in VIG
}


/******************************************************************************\
*			Total number of edges in the VIG								   *
\******************************************************************************/
int Mk::n_edges_VIG( ){
	int n_edges;
	
	for (int i=0;i<N;i++){
		n_edges+=VIG[i].size();	
	}
	
	return (n_edges/2);   		// remember VIG is a simmetric graph
}


/******************************************************************************\
*								Print Mk information						   *
\******************************************************************************/
void Mk::print(void){
	int aux;
	
	cout<< "Mk Landscape: "<<endl<<"  M="<<M<<endl<<"  k= "<< k<<endl<<"  N= "<< N<<endl;

	aux=(int) pow(2.0,k);
	
   cout<<"Variables: subfunctions "<<endl;
	for (int i=0;i<N;i++){
		cout<<" x_"<<i<<": ";
		for(list<int>::iterator j = PSI[i].begin(); j != PSI[i].end(); j++) 
			cout<<"  f_"<<*j;
		cout<<endl;
	}
	
	cout<<"VIG: variables "<<endl;
	for (int i=0;i<N;i++){
		cout<<" x_"<<i<<": ";
		for(list<int>::iterator j = VIG[i].begin(); j != VIG[i].end(); j++) 
			cout<<"  x_"<<*j;
		cout<<endl;
	}
	
	cout<<"Subfuncions: variables "<<endl;
	for (int l=0;l<M;l++){
		cout<<" f_"<<l<<": ";
		for(list<int>::iterator j = PHI[l].begin(); j != PHI[l].end(); j++) 
			cout<<"  x_"<<*j;
		cout<<endl;
	}
	
	cout<<"Subfuncions: contributions "<<endl;
	for (int l=0;l<M;l++){
		cout<<" f_"<<l<<":";
		for (int j=0;j<aux;j++)
			cout<<"  "<<Fl[l][j];
		cout<<endl;
	}
	
}


