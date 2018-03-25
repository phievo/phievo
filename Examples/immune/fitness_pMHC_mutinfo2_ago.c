#define NFUNCTIONS 1 //number of  functions computed by the fitness function. should be at least 1 for the fitness
static double result[NTRIES][NFUNCTIONS]; 			//global matrix containing the all the various fitnesses
static double output_array_raw[NTAU*NLIGANDS*NLIGANDS_SELF];				//global matrix encoding for the response of the network considered (as a vector). 
static double output_array[NLIGANDS][NTAU][NLIGANDS_SELF];				//reshaped version of the above. for example: resp_array[i][j] = 1 (0) means that the ith concentration with jth dissociation time does (not) lead to "immune response".
static double prob_array[NTAU][NPARTITION_OUTPUT];			  		//global matrix containing the "probability" of (tau,response) as computed from the integrated response of a network.
static double output_partition[NPARTITION_OUTPUT];		// contains the sections that are used to partition the output values.


/* The philosophy of the fitness function is to take the total concentration of the Outputs 
for all couples (L_i,t_j). From there we compute a histogram of total outputs by marginalizing over all possible ligand concentrations 
(table prob_array) and compute  mutual information as explained in Lalanne Francois 2013.

output_array_raw is computed below from the last step of integration of the equations in compute_output_array
The mapping between the linear array output_array_raw and output_array is:

e.g. with NTAU = 3, NCONC = 2.

(L0,t0) (L0,t1) (L0,t2) (L1,t0) (L1,t1) (L1,t2)

Inverse mapping:
	n_t = n%NTAU
	n_L = floor(n/NTAU)
	
Forward mapping:
	n = n_L*NTAU + n_t
			
where:		n: 		index in NCELLTOT
			n_t:	index in dissociation time
			n_L: 	index in ligand concentration
			
Note that indices always start at 0.

From output_array_raw, we compute output_array containing all triplets (L_i,t_i,S_i) where S_i is an index corresponding to presence of self or not.
Then binning and marginalization over both ligands and self to compute histograms prob_array is done in compute_prob_array( )

Finally, mutual information is computed in mutual_information() below from prob_array

*/



//initializing all the arrays...
void array_initialization(){
    int i,j,k;
    for(i=0;i<NLIGANDS;i++){
        for(j=0;j<NTAU;j++){
            for(k=0;k<NLIGANDS_SELF;k++){
                output_array[i][j][k] = 0;
            }
        }
    }
    for(i=0;i<(NPARTITION_OUTPUT);i++){
        output_partition[i] = 0;
        for(j=0;j<NTAU;j++){
            prob_array[j][i] = 0;
        }
    }
    for(i=0;i<NCELLTOT;i++){
        output_array_raw[i] = 0;
    }
}

void compute_output_array( double history[SIZE][2][NCELLTOT]){
	int i,j;
	for(i=0;i<NCELLTOT;i++){ // Looping through cells (here effectively concentrations and dissociation times).
		for(j=0;j<NOUTPUT;j++){ //there are two outputs, C_N and D_N
		    output_array_raw[i] += history[trackout[j]][1][i]; 
		}
	}
	/*
	// printing history as a check...
	for(i=0;i<NCELLTOT;i++){
	    printf("Cell %d\n\n",i);
	    
	    for(j=0;j<SIZE;j++){
	        printf("%.4e\n",history[j][1][i]);
	    }
	    printf("Output: %.4e\n",output_array_raw[i]);
	    printf("========================\n\n");
	}
	*/
}

double compute_max( double output_array_raw[NCELLTOT]){
	//simply returns the max of output_array_raw
	int i;
	double max = 0;
	max = output_array_raw[0];
	for(i=1;i<NCELLTOT;i++){
	    if(output_array_raw[i]>max){
	        max = output_array_raw[i];
	    }
	}
	return max;
}

double compute_min( double output_array_raw[NCELLTOT]){
	int i;
	double min = 0;
	min = output_array_raw[0];
	for(i=1;i<NCELLTOT;i++){
	    if(output_array_raw[i]<min){
	        min = output_array_raw[i];
	    }
	}
	return min;
}

void normalize_output_history(){
    int i;
    double max_output;
    max_output = compute_max(output_array_raw);
    for(i=0;i<NCELLTOT;i++){
        output_array_raw[i] = output_array_raw[i]/max_output;
    }
}

/* Function that reshapes the response array.*/
void reshape_output_array(){
	int i;
	int n_lig, n_tau, n_self;
	double buffer;
	for(i=0;i<NCELLTOT;i++){		// Filling the reshaped array according to the map shown above.
	    n_lig = (int)(floor( (i - NTAU*NLIGANDS*floor(i/(NTAU*NLIGANDS)) ) /NTAU ) );
	    n_tau = i%NTAU;
	    n_self = (int)(floor(i/(NTAU*NLIGANDS)));
		output_array[n_lig][n_tau][n_self] = output_array_raw[i];
		//output_array[n_lig][n_tau][n_self] contains the rescaled output concentration with ligand corresponding to n_lig, tau to n_tau, with or without self
	}
	
}

int double_equal(double a, double b){
    return fabs(a - b) < TINY;
}

void my_output_partition(){
    int i;
    double pf;
    pf = POW(10000.0,1.0/(NPARTITION_OUTPUT-1));
    for(i=0;i<(NPARTITION_OUTPUT);i++){
        if(i==0){
             output_partition[i] = 1.0/10000.0;
        }
        else{
            output_partition[i] = 1.0/10000.0*POW(pf,i);
        }
    }
}

/* function returning the "probability" matrix from the "response array".*/
// generate histograms by binning using output_partitions prob[j][i] contains the fraction 
// of initial conditions where total output is in bin corresponding to i, 
// with a tau value corresponding to j


void compute_prob_array( ){
	int i, j, k, l;	// Initialization of counter variables.
	
	
	for(j=0;j<NTAU;j++){	// Initializing the probability array.
		for(i=0;i<(NPARTITION_OUTPUT);i++){
			prob_array[j][i] = 0;
		}
	}
	int count = 0;
	for(j=0;j<NTAU;j++){	// Computing the probability array.
		for(i=0;i<NLIGANDS;i++){
		    for(l=0;l<NLIGANDS_SELF;l++){
		        for(k=0;k<(NPARTITION_OUTPUT-1);k++){
		            		            if( (output_partition[k] <= output_array[i][j][l]) && (output_array[i][j][l] < output_partition[k+1]) ){
		            	
				        prob_array[j][k+1] += 1.0/(NTAU*NLIGANDS*NLIGANDS_SELF);
			            count += 1;
			            //printf("output_array[%d][%d][%d] = %.2e\n",i,j,l,output_array[i][j][l]);
			            //printf("counter1 = %d\n",count);
			        }
			        
			        else if(k==0){
			            if (output_array[i][j][l] < output_partition[0]) {
			                prob_array[j][0] += 1.0/(NTAU*NLIGANDS*NLIGANDS_SELF);
			                count += 1;
			                //printf("output_array[%d][%d][%d] = %.2e\n",i,j,l,output_array[i][j][l]);
			                //printf("counter1 = %d\n",count);
			                //printf("possibility 1\n");
			            }
			        }
			    }
			    if(double_equal(output_partition[NPARTITION_OUTPUT-1], output_array[i][j][l]) ){
				    prob_array[j][NPARTITION_OUTPUT-1] += 1.0/(NTAU*NLIGANDS*NLIGANDS_SELF);
			        count += 1;
			        //printf("output_array[%d][%d][%d] = %.2e\n",i,j,l,output_array[i][j][l]);
			        //printf("counter1 = %d\n",count);
			        //printf("possibility 2\n");    
			    }
			}
		}
	}
	//if(count != NTAU*NLIGANDS*NLIGANDS_SELF){
	//    printf("Number of points binned = %d / %d. Faulty integrator: %d.\n\n",count,NTAU*NLIGANDS*NLIGANDS_SELF,ID_INTEGRATOR);
	//}
	
	/*
	for(k=0;k<NLIGANDS_SELF;k++){
	
	double conc;
	
    printf("\n\n");
	printf("Normalized Output Concentration - L self = %d\n",LIGAND_SELF_LIST[k]);
	printf("-----------------------------------------------------------\n");
	printf("TAU \t\t|");
	for(i=0;i<NTAU;i++){
	    printf("%2.1f \t\t|",TAU_LIST[i]);
	}
	printf("\n");
	printf("-----------------------------------------------------------\n");
	for(i=0;i<NLIGANDS;i++){
	    if(i==0){
	        printf("L \t| %d \t|",LIGAND_LIST[0]);
	    }
	    else{
	        printf(" \t| %d\t|",LIGAND_LIST[i]);
	    }
	    for(j=0;j<NTAU;j++){	        
	        printf("%.2e \t|",output_array[i][j][k]);
	    }
	    printf("\n");
	}
	printf("-----------------------------------------------------------\n");
	
	
	}
	
	printf("\n\n");
	printf("Probability Array\n");
    printf("-------------------------------------------------\n");
    printf("Output (=<) \t|");
    for(i=0;i<NTAU;i++){
        printf("tau_%d = %.1f \t|",i+1,TAU_LIST[i]);
    }
    printf("\n-------------------------------------------------\n");
    for(i=0;i<(NPARTITION_OUTPUT);i++){
        printf("%.2e \t|",output_partition[i]);
        for(j=0;j<NTAU;j++){
            printf("%.4f \t|",prob_array[j][i]);
        }       
        printf("\n");
    }
    printf("-------------------------------------------------\n");
    printf("\n\n");
   
*/       
    
}

/* Computed the mutual information for joint probability distribution given by prob_array (assumed normalized here).
The number of row will be NTAU and the number of column 2 (on or off only).*/
double mutual_information(int ncol, int nrow, double prob_array[][NPARTITION_OUTPUT]){
	
	int i, j;		// Initialization.
	double p_x[nrow], p_y[ncol];
	
	for(i=0;i<nrow;i++){	// Marginal probability for row variable.
		p_x[i] = 0;
		for(j=0;j<ncol;j++){
			p_x[i] += prob_array[i][j];
		}
	}

	for(i=0;i<ncol;i++){	// Marginal probability for column variable.
		p_y[i] = 0;
		for(j=0;j<nrow;j++){
			p_y[i] += prob_array[j][i];
		}
	}
	
	double mut_array[nrow][ncol];
	double mutual_information;
	mutual_information = 0;
	for(i=0;i<nrow;i++){
		for(j=0;j<ncol;j++){
			if(prob_array[i][j]!=0){
				mut_array[i][j] = prob_array[i][j]*log(prob_array[i][j]/p_x[i]/p_y[j]);
			}
			else{
				mut_array[i][j] = 0;
			}
			mutual_information += mut_array[i][j];
		}
	}
	return mutual_information; 
}

static double NOGOODFLAG=0;
void nogood(int ntry){
  //dummy function to fill result if there is a problem (nan, negative concentration)
  NOGOODFLAG=1;
}



/* The fitness function proper. This is the mutual information augmented together with a term preventing detrimental antagonism.*/
void fitness(int k){
    	
    int n,i;
   //We first remove dummy integrations
   result[0][k] = 0;
    for(n=0; n<NCELLTOT; n++)  {
      for (i=0;i<SIZE;i++){
	if (history[i][1][n]<0)  {
	  nogood(k);
	}
	if (history[i][1][n]<1e10)
	  {}
	else {
	  nogood(k);   //intercepts nan
	}
      }
    }
    if (NOGOODFLAG==0)
    	{
    
   		double min_concentration;
    	double max_output_concentration;
    	array_initialization(); //put 0 in all arrays
		compute_output_array(history);
		max_output_concentration = compute_max(output_array_raw);
		my_output_partition(); //generates the logarithmic bins
    	normalize_output_history();//divides by max concentration
		reshape_output_array();//puts output_array_raw into output_array.
		compute_prob_array();// 
	
		min_concentration = 0.1;
    	if(max_output_concentration > min_concentration){
	    	result[0][k] = -(0.95 + rand()%100 * 0.001) *mutual_information(NPARTITION_OUTPUT,NTAU,prob_array)/0.693147180559945; //normalizing by the maximum possible value log(2).
		}
	
	}
}
