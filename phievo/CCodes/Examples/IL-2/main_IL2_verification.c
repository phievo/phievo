

/*  General main.c, should not need any others, see comments in fitness_template for 
the functions that must be supplied in that file for the main to work,
*/

int main()  {
    srand( SEED );
    int i, k,l;
    double score = 0;

    /* Return when no outputs for fitness function */
    if(NOUTPUT <= 0) {
    	printf("%s","no output variables? terminating without integration." );
	return;
    }

    /* Important part of the loop. Integrating the network equations and computing fitness.*/
    for (k=0; k<NTRIES; k++){
    	
    	integrator(k);
	    fitness(k);
	    
	    //need to print something for the rest of the algorithm to rank the networks in the population etc.
	    // CAREFUL: THIS ONLY MAKES SENSE IF NTRIES IS 1!
	    printf("\n\n\n");
	    printf("Doses\n");
	    for(l=0;l<N_DOSES;l++){
	        printf("%.3e\n",DOSES[l]);
	    }
	    printf("\n");
	    printf("N_cell\n");
	    for(l=0;l<N_DOSES;l++){
	        printf("%.3e\n",NUMBER_CELLS[l]);
	    }
	    
	    printf("\n");
	    printf("Norm. Max. Output\n");
	    for(l=0;l<N_DOSES;l++){
	        for(i=0;i<N_N_CELL_;i++){
	            printf("%.3e\t",MAX_OUTPUT_ARRAY[l][i]);
	        }
	        printf("\n");
	    }
	    printf("\n");
	    
	    // For better readability.
	    // Printing the max_output and all that to see if everything makes sense.
	    printf("\n");
	    printf("Norm. Max. Output\n");
	    printf("dose\\N\t|\t");
	    for(i=0;i<N_N_CELL_;i++){
	        printf("%.1e\t",NUMBER_CELLS[i]);
	    }
	    printf("\n");
	    printf("---------------");
	    for(i=0;i<N_N_CELL_;i++){
	        printf("--------");
	    }
	    printf("\n");
	    for(l=0;l<N_DOSES;l++){
	        printf("%.1e\t|\t",DOSES[l]);
	        for(i=0;i<N_N_CELL_;i++){
	            printf("%.3f\t",MAX_OUTPUT_ARRAY[l][i]);
	        }
	        printf("\n");
	    }
	    printf("\n");
	    
	    
	    
	    // Printing the probability array for the input
	    printf("\n");
	    printf("Prob. Array. Dose\n");
	    printf("dose\\Output|\t");
	    for(i=0;i<N_PARTITION_OUTPUT;i++){
	        printf("%.1e\t",output_partition[i]);
	    }
	    printf("\n");
	    printf("---------------");
	    for(i=0;i<N_PARTITION_OUTPUT;i++){
	        printf("--------");
	    }
	    printf("\n");
	    for(l=0;l<N_DOSES;l++){
	        printf("%.1e\t|\t",DOSES[l]);
	        for(i=0;i<N_PARTITION_OUTPUT;i++){
	            printf("%.3f\t",prob_array_dose[l][i]);
	        }
	        printf("\n");
	    }
	    printf("I(out;in) = %.3f\n",-result[0][0]);
	    
	    
	    // Printing the probability array for the input
	    printf("\n");
	    printf("Prob. Array. N_cell\n");
	    printf("N_cell\\Output|\t");
	    for(i=0;i<N_PARTITION_OUTPUT;i++){
	        printf("%.1e\t",output_partition[i]);
	    }
	    printf("\n");
	    printf("---------------");
	    for(i=0;i<N_PARTITION_OUTPUT;i++){
	        printf("--------");
	    }
	    printf("\n");
	    for(l=0;l<N_N_CELL_;l++){
	        printf("%.1e\t|\t",NUMBER_CELLS[l]);
	        for(i=0;i<N_PARTITION_OUTPUT;i++){
	            printf("%.3f\t",prob_array_ncell[l][i]);
	        }
	        printf("\n");
	    }
	    printf("I(out;N_cell) = %.3f\n\n",result[0][1]);
	    
	    
       	if( PRINT_BUF ){
       	    for(l=0;l<NCELLTOT;l++){
	            print_history(k,l);
	        }
	    }
    }
}  

