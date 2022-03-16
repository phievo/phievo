/*
set history of all species at initial time to rand number in [0.1,1) eg avoid 0 and initializes Input.
For this problem, we created a function init_signal that just make a random gate function.

tfitness is defined in the fitness C file
*/
static double isignal[NSTEP][NCELLTOT][2];
 int next_time(){
        return 100+(rand()%500);
    }
void init_signal( ){
   
    
    int k, t, tnext, val;
    // Construct the first signal
    for(k=0; k<NCELLTOT; k++){
        tnext = next_time();
        val = rand()%2;
        for(t=0; t<NSTEP; t++){
            tnext -= 1;
            isignal[t][k][0] = val;
            if (tnext <= 0){
                tnext = next_time();
                val = (val+1)%2;
            }
        }
    }
    // Construct the second signal
    for(k=0; k<NCELLTOT; k++){
        tnext = next_time();
        val = rand()%2;
        for(t=0; t<NSTEP; t++){
            tnext -= 1;
            isignal[t][k][1] = val;
            if (tnext <= 0){
                tnext = next_time();
                val = (val+1)%2;
            }
        }
    }
}

void init_history(int kk){
    init_signal();
    int ncell,n_gene;
    for (ncell=0;ncell<NCELLTOT;ncell++){
        for (n_gene=0;n_gene<SIZE;n_gene++){
            history[n_gene][0][ncell] = 0;
        }
    }
}
