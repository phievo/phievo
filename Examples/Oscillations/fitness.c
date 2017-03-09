

#define NFUNCTIONS 1 //number of  functions computed by the fitness function. should be at least 1 for the fitness
#define NGENE 3
static double result[NTRIES][NFUNCTIONS]; //result will contain the fitness plus the other results computed by the fitness function

static int tfitness=NSTEP/6;//time from which we start pulses of I and compute the fitness

void nogood(int ntry){
  //dummy function to return infinite results if some fitness criteria - like crazy concentrations - is realized
  result[ntry][0]=RAND_MAX;
}

void fitness(double history[][NSTEP][NCELLTOT], int trackout[],int ntry){
  int t,i, g,  prod, bin;
  int input = trackin[0];
  
  double fitness = 0;
  double count_accept  = 0;
  double probG[NSTEP][NGENE];
  double probI[NGENE];
  double bining = 1.0;
  double probIG[NGENE][NGENE];
  

  
  double sumI=0,sumG[NSTEP],sumIG=0;
  for(t=0;t<NSTEP;t++)sumG[t]=0;
  for(t=0;t<NSTEP;t++){
    
    int indexI = (int)(history[input][t][0]/bining);
    
    probI[indexI] += 1;
    
    sumI+=1;
    for(g=0;g<NGENE;g++){
      probG[t][g] += history[trackout[g]][t][0];
      sumG[t] += probG[t][g];
      probIG[indexI][g] += history[trackout[g]][t][0];
      sumIG += probIG[indexI][g];
    }        
  }
  
  for(g=0;g<NGENE;g++){
    //probI[g] /= sumI;    
    for(i=0;i<NGENE;i++){
      probIG[i][g] /= sumIG;
      fitness += probIG[i][g]*log(probIG[i][g]);
    }
    for(t=0;t<NSTEP;t++){
      probG[t][g] /= sumG[t];
      fitness -= probG[t][g]*log(probG[t][g]);
    }
  }
  if(fitness != fitness){
    fitness = RAND_MAX;
  }
  result[ntry][0] = fitness;
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
