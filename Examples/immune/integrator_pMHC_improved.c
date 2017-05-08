

 /* Takes y_n to y_{n+step} according to RK4. */
void RK(double y_n[SIZE], double step, double tau){
    double k1[SIZE], k2[SIZE], k3[SIZE], k4[SIZE], temp[SIZE];

    int i;
	for(i=0;i<SIZE;i++){
    	k1[i]=0.0;
    	k2[i]=0.0;
    	k3[i]=0.0;
    	k4[i]=0.0;
    	temp[i]=0.0;
    }

    derivC(y_n,k1,tau,step);			// Modifies k1 to hold the first approximation of d(y_n) and outputs the noise_increment.
    for(i=0;i<SIZE;i++){		// temp is the input concentrations for the next computation.
        temp[i] = y_n[i] + 0.5*step*k1[i];
    }

    derivC(temp,k2,tau,step);	// Modifies k2
    for(i=0;i<SIZE;i++){
        temp[i] = y_n[i] + 0.5*step*k2[i];
    }

    derivC(temp,k3,tau,step); 	// Modifies k3.
    for(i=0;i<SIZE;i++){
        temp[i] = y_n[i] + step*k3[i];
    }

    derivC(temp,k4,tau,step);	// Modifies k4.

    for(i=0;i<SIZE;i++){	// Puts it all together.
        y_n[i] += step*1.0/3.0*(0.5*k1[i] + k2[i] + k3[i] + 0.5*k4[i]);
    }
}

/* Takes the tolerance, scale and difference vector and returns the maximum ratio of (difference)/(desired tolerance) */
double desired_tolerance(double tol, double y_scale[SIZE], double diff[SIZE], int *ptrID){
    int i;
    double max_ratio;
    max_ratio = diff[0]/tol/y_scale[0]; // We wish to store the maximum of this ratio to evaluate the next step to take in the algorithm.
    *ptrID = 0;
    for(i=1;i<SIZE;i++){
    	if( max_ratio < (diff[i]/tol/y_scale[i]) ){
    	    max_ratio = diff[i]/tol/y_scale[i];		// If the so far maximum ratio is smaller than a value, actualize the maximum ratio to this new value.
    	    *ptrID=i;
    	}
    }
    return max_ratio;
}

/* Moves forward using one and two steps and compare the results to adapt the time step.
   Returns the next time step and actualizes the value of y_before to y_after...
   A pointer to a time scalar and number of succeeded step must be provided.*/
double RKstepper(double y_before[SIZE], double y_scale[SIZE], double tol, double *ptrstep, double tau, int *ptrnstep, int ncell){

    // Initializing arrays containing the temporary solutions of the local integration.
    double temp1[SIZE], temp2[SIZE];

    int i;
    for(i=0;i<SIZE;i++){
        temp1[i] = y_before[i];
        temp2[i] = y_before[i];
    }

    // Evaluating the first step in one shot. Storing the result in temp1.
    RK(temp1,*ptrstep,tau);

    // Evalutating the same step but in two smaller steps of step/2.
    RK(temp2,0.5*(*ptrstep),tau);
    RK(temp2,0.5*(*ptrstep),tau);

    // Evaluating the difference between the two to assess the error in our procedure.
    double diff[SIZE];
    for(i=0;i<SIZE;i++){
        diff[i] = fabs(temp2[i] - temp1[i]);
    }

    // Next step depends if the required tolerance is met.
    double max_ratio;
    int ID;
    max_ratio = desired_tolerance(tol,y_scale,diff, &ID);

    // If the tolerance condition is met.
    if (max_ratio<1){

    	// Actualizing the value of the concentrations (with noise)
    	for(i=0;i<SIZE;i++){
    	    y_before[i] = temp2[i];
    	}

    	// Actualizing the next time step to try in the algorithm (the increase in step size is bounded by ABOVE).
    	if(POW(max_ratio,-0.2)<ABOVE){
    	    *ptrstep = (*ptrstep)*POW(max_ratio,-0.2);
    	}
    	else{
    	    *ptrstep = (*ptrstep)*ABOVE;
    	}

    }
    // If the tolerance condition is not met, calls recursively RKstepper with the smaller step.
    else{

		// Estimating the next appropriatee time step (the decrease in step size is bounded by BELOW).
        if(POW(max_ratio,-0.25) > BELOW){
            *ptrstep = 0.7*(*ptrstep)*POW(max_ratio,-0.25);
        }
        else{
            *ptrstep = (*ptrstep)*BELOW;
        }

        RKstepper(y_before,y_scale,tol,ptrstep,tau,ptrnstep,ncell);
    }

    	// moving forward in time
    	return *ptrstep;
}



void myODEINT(double max_time, double history[SIZE][2][NCELLTOT], double init_step, int ncell, double tol){
    int nstep = 0;
    int i;
    double output = 0;
    double integrated_species[SIZE];

    // The vector storing the present state of the system. Initialized to conditions at t=0.
    // Initializing the integrated concentration array.
    double s[SIZE];
    for(i=0;i<SIZE;i++){
         s[i] = history[i][0][ncell];
         integrated_species[i] = 0;
    }

    // We take an absolute error scale for the very first step.
    double scale[SIZE];
    for(i=0;i<SIZE;i++){
        scale[i] = 1;
    }

    double time = 0;  // Initializing the time variable at each new cell (i.e. for each new initial conditions).
    double step = init_step;

    while(time < max_time){

		// Takes a step forward in time with trial step size "step" (which is modified dynamically by the algorithm).
        step = RKstepper(s,scale,tol,&step,N_to_TAU(ncell),&nstep,ncell);
        time += step;
        for(i=0;i<SIZE;i++){
        	scale[i] = s[i] + TINY;				// we take the scale to be relative and avoid division by 0 by adding a small number.
        }

        // Checking for negative concentrations and integrating the output.
        for(i=0;i<SIZE;i++){
            if(s[i] < 0){   //negative concentration is possible in Langevin (we do the crude thing of putting such negative concentrations to 0...).
                s[i] = 0;
            }
            integrated_species[i] += step*s[i];

        }
        //output += step*s[trackout[0]];
    }
    // Filling the history array with only the final steady-state values.
    for(i=0;i<SIZE;i++){
        history[i][1][ncell] = integrated_species[i]/time;

        /*
        if(i == trackout[0]){
            history[i][1][ncell] = output/time;
        }
        else{
            history[i][1][ncell] = s[i];
        }
        */
    }

}


void integrator(int kk){

    int ncell;
    init_history(kk);

    /* Integrates the network equations for each "cell" using adaptive Runge-Kutta of order 4. */
    for(ncell=0;ncell<NCELLTOT;ncell++){
        //printf("Now at ncell = %d.\n",ncell);
        myODEINT(MAX_TIME,history,DT,ncell,TOL);
    }
}
