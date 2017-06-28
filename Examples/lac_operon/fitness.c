 /*
Defines the fitness function for a logical gate fitness
Last Edited the 06 feb. 2017
Coder : M. Hemery
*/

#define NFUNCTIONS 1 //number of  functions computed by the fitness function. should be at least 1 for the fitness

static double result[NTRIES][NFUNCTIONS]; //result will contain the fitness plus the other results computed by the fitness function

static int tfitness=NSTEP/6;//time from which we start pulses of I and compute the fitness

void nogood(int ntry){
    //dummy function to return infinite results if some fitness criteria - like crazy concentrations - is realized
    result[ntry][0]=RAND_MAX;
}

void fitness(double history[][NSTEP][NCELLTOT], int trackout[],int ntry){
    int t, prod;
    double conc;
    int input0 = trackin[0];
    int input1 = trackin[1];
    int output = trackout[0];
    double score = 0;
    double best = 0;
    
    // Collect direct data
    for (t=0; t<NSTEP; t++){
        prod = (int) (history[input0][t][0] * history[input1][t][0]);
        conc = history[output][t][0];
        if(conc>1.){conc = 1.;}
        score += (1.75*prod - .75)*conc;
        best += prod;
    }
    
    result[ntry][0] = score/best;
}

void treatment_fitness( double history2[][NSTEP][NCELLTOT], int trackout[]){
    /* function to print out the result*/
    //if you want to do anything with the average output history2, this is the right place
    int k;
    double mean;
    // Compute the mean
    mean = 0;
    for (k=0; k<NTRIES; k++){
        mean += result[k][0];
    }
    mean /= NTRIES;
    printf("%f",-mean);
}
