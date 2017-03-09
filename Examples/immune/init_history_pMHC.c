
/* set history of all species at initial time to rand number in [0.1,1) eg avoid 0 */
int contains(int array[pMHC_LENGTH],int element,int size){
    int i;
    for(i=0;i<size;i++){
        if(array[i]==element){
            return 1;
        }
    }
    return 0;
}

/*function returning the ligand concentration index corresponding to a given cell number. */
int N_to_L(int ncell)  {
    declare_list();
    //if(ncell<(NTAU*NLIGANDS)){
        return LIGAND_LIST[(int)(floor( (ncell - NTAU*NLIGANDS*floor(ncell/(NTAU*NLIGANDS)) ) /NTAU ) )];
    //}
    //else{
    //    return LIGAND_LIST[ncell-(NTAU*NLIGANDS)];
    //}
}

/*function returning the ligand concentration corresponding to a given cell number. */
double N_to_TAU(int ncell)  {
    declare_list();
    //if(ncell<(NTAU*NLIGANDS)){
        return TAU_LIST[ncell%NTAU];
    //}
    //else{
    //    return TAU_LIST[NTAU-1];		// returning the largest value of binding time for the agonist time.
    //}
}


double N_to_L_SELF(int ncell)  {
        return LIGAND_SELF_LIST[(int)(floor(ncell/(NTAU*NLIGANDS)))];
}



void init_history()  {
    int ncell,n_species;
    for (ncell=0;ncell<NCELLTOT;ncell++){
    	for (n_species=0;n_species<SIZE;n_species++){
    	    if(n_species==0){
      	        history[n_species][0][ncell] = N_to_L(ncell); 	// this is the initial free ligand concentration.
      	    }
      	    else if(n_species==1){
      	        history[n_species][0][ncell] = RECEPTOR;		// this is the initial free receptor concentration.
      	    }
      	    else if(contains(pMHC_LIST,n_species,pMHC_LENGTH)){
      	        history[n_species][0][ncell] = 0;   			// elements of the pMHC cascade initialized to 0.
      	    }
      	    else if(contains(KP_LIST,n_species,KP_LENGTH)){		// concentration of kinases or phosphatase (fixes it between generations and allows to put selection pressure on the concentration).
      	        history[n_species][0][ncell] = KP_CONC[n_species];
      	    }
      	    else if(n_species==(SIZE-pMHC_LENGTH-1)){			// concentration of self ligands present
      	        history[n_species][0][ncell] = N_to_L_SELF(ncell);
      	    }
      	    //else if(n_species>SIZE-pMHC_LENGTH-1){
      	    //    history[n_species][0][ncell] = 0;
      	    //}
      	    else{												// this takes into account phosphorylated kinase/phosphatase as well as pMHC complexes associated with antagonistic ligands.
      	        history[n_species][0][ncell] = 0;
      	    }
      	}
    }
}

