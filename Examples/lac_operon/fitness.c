 /*
Defines the fitness function for the exponential evolution
Create the 10 march 2016 from the Adaptation fitness
Last Edited the 21 march 2016
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
    int t, prod, bin;
    int input0 = trackin[0];
    int input1 = trackin[1];
    int output = trackout[0];
    
    double data[2][20];
    double dataprod[2];
    double databin[20];
    
    for (bin=0; bin<20; bin++){
        data[0][bin] = 0;
        data[1][bin] = 0;
        databin[bin] = 0;
    }
    dataprod[0] = 0;
    dataprod[1] = 0;
    
    // Collect direct data
    for (t=0; t<NSTEP; t++){
        prod = ((int) history[input0][t][0]+1)%2 * history[input1][t][0];
        bin = (int) history[output][t][0]*10;
        if(bin>=20){bin = 19;}
        data[prod][bin] += 1;
        dataprod[prod] += 1;
        databin[bin] += 1;
    }
    
    // Compute mutual information
    double MI = 0;
    double Ppb,Pb,Pp;
    for (bin=0; bin<20; bin++){
        for (prod=0; prod<2; prod++){
            if (data[prod][bin] > 0){
                Ppb = (double) data[prod][bin]/(double)NSTEP;
                Pp = (double) dataprod[prod]/(double)NSTEP;
                Pb = (double) databin[bin]/(double)NSTEP;
                MI += (double)  Ppb*log(Ppb/(Pp*Pb));
            }
        }
    }
    result[ntry][0] = MI;
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
    mean /= NTRIES*.5623;
    printf("%f",-mean);
}
