// Store the results for individual trials.

static double result[NTRIES];

void fitness( double history[][NSTEP][NCELLTOT], int trackout[],int ntry)  {
  int trial = 0;
  int func = 0;
  for(trial=0;trial<NTRIES;trial++)
    {
      result[trial] = 0;
	
    }
}

// Combine the fitnesses obtained in the different trials.
void treatment_fitness( double history2[][NSTEP][NCELLTOT], int trackout[]){
  int trial = 0;
  int func = 0;
  double total_fitness = 0;
  for(trial=0;trial<NTRIES;trial++)
    {
      total_fitness += result[trial];
    }
  
  // Print network's fitness caught by the python code
  printf("%f",total_fitness);
}
