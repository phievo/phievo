/* set history of all species at initial time to rand number in [0.1,1) eg avoid 0 and initializes Input.
For this problem, we created a function init_signal that defines random steps of inputs.
init_signal creates time dependent input from [0, NSTEP) by integrating a sum of random positive delta functions and then exponentiating.

tfitness is defined in the fitness C file
*/

static double isignal[NSTEP][NCELLTOT];
static double dsignal[NSTEP][NCELLTOT];


static int t1 = 1;   // duration of delta function spike in dsignal, dt units
static int t2 = 1000;  // mean time between delta fns (uniform distrib [1, 2*t2]  dt units


void init_signal( ) {

    int k, t, tnext;
    double ampl;	// random value of dsignal

    for(k=0; k<NCELLTOT; k++)  {
      for(t=0; t<NSTEP; t++) {
	 dsignal[t][k] = 0;   //if multiple calls need zero expl
      }
    }


    tnext=tfitness;

    for(k=0; k<NCELLTOT; k++)  {

      while( tnext < NSTEP ) {
	tnext += 500+(int) 2*t2*(rand()/( (double) RAND_MAX +1));
	ampl = 20*(2*(rand()/((double)RAND_MAX + 1)) -1);

	if(tnext + t1 > NSTEP) {
	  break;
	}
	for( t=tnext; t<tnext+t1; t++) {
	  dsignal[t][k] = ampl;
	}
	tnext += t1;
      }
    }

/* integrate the derivative and expon to get real signal applied to cell  */

    for(k=0; k<NCELLTOT; k++)  {
      isignal[0][k] = (2*(rand()/((double)RAND_MAX + 1)) -1);
      for(t=1; t<NSTEP; t++) {
	isignal[t][k] = isignal[t-1][k] + DT*dsignal[t][k];
      }
      for(t=0; t<NSTEP; t++) {
	isignal[t][k] = exp(isignal[t][k]);
      }
    }

}


void init_history(int trial)  {

    int ncell,n_gene;
    init_signal();  // initialize signal for all times.

    for (ncell=0;ncell<NCELLTOT;ncell++){
    	for (n_gene=0;n_gene<SIZE;n_gene++){
	  history[n_gene][0][ncell] = FRAND();
        }
    }
}
