/******************************************************************************\
*								Perturbation								   *
\******************************************************************************/
#include "defs.h"


/******************************************************************************\
*	Perturbation based on randomly kicking rball_size bits			  		   *
\******************************************************************************/
void randomKickbP(int *x, int *x_new, int rball_size){
	int r, *order_v; 
	
	for (int i=0;i<chrom_size;i++)
		x_new[i]=x[i];
	
	if (rball_size==0 || rball_size>chrom_size)
		return;

	// Randomly defining rball_size bits 
	order_v=aloc_vectori(chrom_size);	
	rand_perm_size(vsort_aux, order_v, chrom_size, rball_size);
	
	// Randomly changing rball_size bits		
	for (int i=0;i<rball_size;i++){
		r=order_v[i];
		
		if (x_new[r]==0)
			x_new[r]=1;
		else
			x_new[r]=0;
	}


	delete [] order_v;
	
	return;
			
}	
