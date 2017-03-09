
/* set history of all species at initial time to rand number in [0.1,1) eg avoid 0 */
int contains(int array[],int element,int size){
    int i;
    for(i=0;i<size;i++){
        if(array[i]==element){
            return 1;
        }
    }
    return 0;
}

/*function returning the number of cells as function of trial label (ncell). */
int N_to_N_cell_(int ncell)  {
    declare_list();
    return NUMBER_CELLS[(int)(floor(ncell/N_DOSES))];
}

/*function returning the input dose as function of trial label (ncell). */
double N_to_Dose(int ncell)  {
    declare_list();
    return DOSES[ncell%N_DOSES];
}


// Note that the concentration of the input will be overwritten in the integrator. 
void init_history()  {
    int ncell,n_species;
    for (ncell=0;ncell<NCELLTOT;ncell++){
    	    	
    	for (n_species=0;n_species<SIZE;n_species++){
      	    if(contains(IC_LIST,n_species,IC_LENGTH)){		// concentration of kinases or phosphatase (fixes it between generations and allows to put selection pressure on the concentration).
      	        history[n_species][0][ncell] = IC_CONC[n_species];
      	    }
      	    else{												
      	        history[n_species][0][ncell] = 0;
      	    }
      	}
    }
}

