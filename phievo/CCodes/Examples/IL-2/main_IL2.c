

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
	return 0;
    }

    /* Important part of the loop. Integrating the network equations and computing fitness.*/
    for (k=0; k<NTRIES; k++){
    	
    	integrator(k);
	    fitness(k);
	    
	    //need to print something for the rest of the algorithm to rank the networks in the population etc.
	    // CAREFUL: THIS ONLY MAKES SENSE IF NTRIES IS 1!
	    printf("%f\n",result[0][0]);
	    printf("%f",result[0][1]);
	    
       	/*if( PRINT_BUF ){
	        print_history(k);
	    }*/
    }
return 1;
}  

