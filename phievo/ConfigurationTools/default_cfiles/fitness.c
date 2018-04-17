// Fitness function
// !!! THIS FILE MUST BE UPDATED WITH AN ACTUAL FITNESS FUNCTION !!!
// For more details:
// http://phievo.readthedocs.io/en/latest/create_new_project.html#fitness-c
// Examples of fitnesses can be found in the example projects:
// https://github.com/phievo/phievo/tree/master/Examples
// Store the results for individual trials.
static double result[NTRIES];

void fitness( double history[][NSTEP][NCELLTOT], int trackout[],int trial)  {
  // MODIFY HERE to compute a fitness from the history array
  result[trial] = 1; 
}

// Combine the fitnesses obtained in the different trials.
void treatment_fitness( double history2[][NSTEP][NCELLTOT], int trackout[]){
  // MODIFY HERE to combine the fitnesses computed for the different trials.
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
