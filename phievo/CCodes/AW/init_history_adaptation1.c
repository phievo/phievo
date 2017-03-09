
/* set history of all species at initial time to rand number in [0.1,1) eg avoid 0 */

void init_history()  {

    int ncell,n_gene;
    init_signal();  // initialize signal for all times.

    for (ncell=0;ncell<NCELLTOT;ncell++){
    	for (n_gene=0;n_gene<SIZE;n_gene++){
	  history[n_gene][0][ncell] = FRAND();
        }
    }
}

