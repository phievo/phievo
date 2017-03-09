static double result[NTRIES][NFUNCTIONS]; 						//global matrix containing the all the various fitnesses
static double MAX_OUTPUT_ARRAY[N_DOSES][N_N_CELL_];	//reshaped version of the MAX_OUTPUT array.
static double prob_array_dose[N_DOSES][N_PARTITION_OUTPUT];			  	//global matrix containing the "probability" of (tau,response) as computed from the integrated response of a network.
static double prob_array_ncell[N_N_CELL_][N_PARTITION_OUTPUT];
static double output_partition[N_PARTITION_OUTPUT];				// contains the sections that are used to partition the output values.

/* The fitness function takes the max concentration of the Outputs 
for all couples (n_d_j,n_c_k), where n_d_j is a total dose of ligand and n_c_k number of cells exposed to it.
 This is first stored in linear array MAX_OUTPUT, computed when we integrate equations.
We then map linear array MAX_OUTPUT to a matrix MAX_OUTPUT_ARRAY so that MAX_OUTPUT_ARRAY[n_d_j,n_c_k] contains corresponding output concentration
The mapping between the linear array MAX_OUTPUT  and MAX_OUTPUT_ARRAY is:

e.g. with  N_N_CELL = 2 N_DOSES = 3

MAX_OUTPUT[i]
0:(n_d_0,n_c_0) 1:(n_d_1,n_c_0) 2:(n_d_2,n_c_0) 3:(n_d_0,n_c_1) 4:(n_d_1,n_c_1) (n_d_2,n_c_1) 

Inverse mapping:
	
		n_dose = i%N_DOSES;
	    n_cell_ = (int)(floor(i/N_DOSES));
Forward mapping:
	n = n_cell_*N_DOSES + n_dose
			
where:		n: 		index in NCELLTOT
			n_t:	index in dissociation time
			n_L: 	index in ligand concentration
			
Note that indices always start at 0.

Then binning and marginalization over doses/cells to compute histograms prob_array_dose and prob_array_ncell is done in compute_prob_arrays( ). 
Note the binning of output is done on log scale.
prob_array_dose[N_DOSES][N_PARTITION_OUTPUT]: histograms of Max output marginalizing over cell number, i.e. when only doses vary
prob_array_dose[N_N_CELL_][N_PARTITION_OUTPUT]: histograms of Max output marginalizing over doses, i.e. when only cell numbers vary
Then we compute mutual information  to see how much information the Max output contains respectively on the ligand concentration and on the number of cells
in mutual_information() from prob_arrays


*/



//initializing all the arrays...
void array_initialization(){
    int i,j,k;
    for(i=0;i<(N_PARTITION_OUTPUT);i++){
        output_partition[i] = 0;
        for(j=0;j<N_DOSES;j++){
            prob_array_dose[j][i] = 0;
        }
        for(k=0;k<N_N_CELL_;k++){
            prob_array_ncell[k][i] = 0;
        }
    }
}

double compute_max( double MAX_OUTPUT[NCELLTOT]){
	int i;
	double max = 0;
	max = MAX_OUTPUT[0];
	for(i=1;i<NCELLTOT;i++){
	    if(MAX_OUTPUT[i]>max){
	        max = MAX_OUTPUT[i];
	    }
	}
	return max;
}

double compute_min( double MAX_OUTPUT[NCELLTOT]){
	int i;
	double min = 0;
	min = MAX_OUTPUT[0];
	for(i=1;i<NCELLTOT;i++){
	    if(MAX_OUTPUT[i]<min){
	        min = MAX_OUTPUT[i];
	    }
	}
	return min;
}

void normalize_output_history(){
	//divide by max output to have all MAX_OUTPUT below 1
    int i;
    double max_output;
    max_output = compute_max(MAX_OUTPUT);
    for(i=0;i<NCELLTOT;i++){
        MAX_OUTPUT[i] = MAX_OUTPUT[i]/(max_output + TINY);
    }
}

/* Function that reshapes the response array.*/
void reshape_output_array(){
	int i;
	int n_dose, n_cell_;
	double buffer;
	for(i=0;i<NCELLTOT;i++){		// Filling the reshaped array according to the map shown above.
	    n_dose = i%N_DOSES;
	    n_cell_ = (int)(floor(i/N_DOSES));
		MAX_OUTPUT_ARRAY[n_dose][n_cell_] = MAX_OUTPUT[i];
	}
	
}

int double_equal(double a, double b){
    return fabs(a - b) < TINY;
}

void my_output_partition(){
	//defines logarithmic bin
    int i;
    double pf;
    pf = POW(10000.0,1.0/(N_PARTITION_OUTPUT-1));
    for(i=0;i<(N_PARTITION_OUTPUT);i++){
        if(i==0){
             output_partition[i] = 1.0/10000.0;
        }
        else{
            output_partition[i] = 1.0/10000.0*POW(pf,i);
        }
    }
}

/* function returning the "probability" matrix from the "response array".*/
void compute_prob_arrays( ){
	int i, j, k, l;	// Initialization of counter variables.
	
	// filling the probability array for the input dose (used to compute I(input;output))
	for(j=0;j<N_DOSES;j++){	// Computing the probability array.
		for(i=0;i<N_N_CELL_;i++){
		    for(k=0;k<N_PARTITION_OUTPUT;k++){
		        
		        if(k==0){
			        if (MAX_OUTPUT_ARRAY[j][i] <= output_partition[k]) {
			            prob_array_dose[j][k] += 1.0/(N_DOSES*N_N_CELL_);
			        }
			    }
		            
		        else if( (output_partition[k-1] < MAX_OUTPUT_ARRAY[j][i]) && ((MAX_OUTPUT_ARRAY[j][i] <= output_partition[k]) || (double_equal(MAX_OUTPUT_ARRAY[j][i],output_partition[k])) )){
				    prob_array_dose[j][k] += 1.0/(N_DOSES*N_N_CELL_);
			    }
			}
		}
	}
	
	// filling the probability array for the number of cells (used to compute I(N_cel;output))
	for(j=0;j<N_N_CELL_;j++){	// Computing the probability array.
		for(i=0;i<N_DOSES;i++){
		    for(k=0;k<N_PARTITION_OUTPUT;k++){
		        
		        if(k==0){
			        if (MAX_OUTPUT_ARRAY[i][j] < output_partition[k]) {
			            prob_array_ncell[j][k] += 1.0/(N_DOSES*N_N_CELL_);
			        }
			    }
			    
		        else if( (output_partition[k-1] < MAX_OUTPUT_ARRAY[i][j]) && ((MAX_OUTPUT_ARRAY[i][j] <= output_partition[k]) || (double_equal(MAX_OUTPUT_ARRAY[i][j],output_partition[k])) )){
				    prob_array_ncell[j][k] += 1.0/(N_DOSES*N_N_CELL_);
			    }
			}
		}
	}
}
	
    

/* Computed the mutual information for joint probability distribution given by prob_array (assumed normalized here).
The number of row will be NTAU and the number of column the number of elements in the output_partition.*/
double mutual_information(int ncol, int nrow, double prob_array[][N_PARTITION_OUTPUT]){
	
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

/* The fitness function proper. This is the mutual information augmented together with a term preventing detrimental antagonism.*/
void fitness(k){
    
    int i;
    double min_concentration;
    double max_output_concentration;
    array_initialization();
    
    // various functions to reshape output arrays, get the probability arrays etc.
	my_output_partition();
    normalize_output_history();
    
    
	reshape_output_array();
	compute_prob_arrays();
	
	// evaluating the two (normalized) mutual informations.
	//result[k][0] = -1.0*mutual_information(N_PARTITION_OUTPUT,N_DOSES,prob_array_1)/log(N_DOSES)/(1.0 +mutual_information(N_PARTITION_OUTPUT,N_N_CELL_,prob_array_2)/log(N_N_CELL_));
	result[k][0] = -(0.95 + rand()%100 * 0.001)*mutual_information(N_PARTITION_OUTPUT,N_DOSES,prob_array_dose)/log(N_DOSES);
	result[k][1] = (0.95 + rand()%100 * 0.001)*mutual_information(N_PARTITION_OUTPUT,N_N_CELL_,prob_array_ncell)/log(N_N_CELL_);	
	
	// NOTE: the above function that is with a -1 is MAXIMIZED, the other MINIMIZED.
	}
	
/*
	result[0][k] = -1.0*mutual_information(N_PARTITION_OUTPUT,N_DOSES,prob_array_1)/log(N_DOSES)/(1.0 +mutual_information(N_PARTITION_OUTPUT,N_N_CELL_,prob_array_2)/log(N_N_CELL_));	// we want to maximize this, this is why we need the -1 in front (the algorithm selects for the minimum...)
*/
