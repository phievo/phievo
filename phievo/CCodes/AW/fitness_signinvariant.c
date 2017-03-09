/* Template for the fitness, the fitness() and treatment_fitness() are used by the
main.c and must be supplied with the arguments and returns indicated.
   The fitness() called for each 'try' (eg iteration over initial conditions and inputs)
Any results needed for computing fitness for selection and other diagonstic to be
returned to python code are accumulated in result[][] array.   
   treatment_fitness is called by main after all ntries are run and 
computes what is to be returned to python.  At this stage can compute the fitness
from history[][][] arveraged over all tries.
*/

#define NFUNCTIONS 1 //number of  functions computed by the fitness function. should be at least 1 for the fitness

static double result[NTRIES][NFUNCTIONS]; //global matrix containing the all the various fitnesses
static long idum = -81091;
    

void nogood(int ntry){
  //dummy function to fill result if there is a problem (nan, negative concentration)

  int i=0;
  for (i=0;i<NFUNCTIONS;i++){
     result[ntry][i]=RAND_MAX;
  }
}
    
void  fitness( double history[][NSTEP][NCELLTOT], int trackout[],int ntry)  {

    
    int k,i;
    for(k=0; k<NCELLTOT; k++)  {
      for (i=0;i<SIZE;i++){
	if (history[i][NSTEP-1][k]<0)  {
	  nogood(ntry);
	}
	if (history[i][NSTEP-1][k]<1e10)
	  {}
	else {
	  nogood(ntry);   //intercepts nan
	}
      }
    }
    /* Compute the fitness here */
	
	int ind0 = trackout[0];
	double midpoint = history[ind0][NSTEP/2][0];
	double x,xa, xmax = 0,xamax = 0;
	for( i = NSTEP/2+1; i < NSTEP; i++){
		x = history[ind0][i][0]-midpoint;
		xa = fabs(x);
		if(xa > xamax){ xamax = xa; xmax = x;}
	}
	result[ntry][0] = xmax;
	
	 

}


void treatment_fitness( double history2[][NSTEP][NCELLTOT], int trackout[]){
 
  /* function to print out the result*/



  //if you want to do anything with the average output history2, this is the right place


  int l,k;
  double score,response,difference;
  double stat_score[NTRIES];
  for (k=0;k<NTRIES;k++){
    stat_score[k]=result[k][0];
  }

  /* define score that will be used for selection here.  For most problems 
  score and average are the same  */
  
  response = fabs(result[0][0])+fabs(result[1][0]);
  difference = fabs(result[0][0]-result[1][0]);
 
  response = -response/(5+response);
 
  score = difference + response;
  

  /*score=average_score(stat_score);
  double average=average_score(stat_score);
  double variance_score=std_dev(stat_score); */

  if (score < RAND_MAX){
      //printf("%f\n%f\n%f\n",score,average,variance_score);     
      printf("%f\n%f\n%f",difference,response,score);  
  }
  else
      printf("%i\n%i\n%i",RAND_MAX,RAND_MAX,RAND_MAX);
  if( NFUNCTIONS > 1) printf("\n");   

  for (l=1;l<NFUNCTIONS;l++){
    for (k=0; k<NTRIES; k++){
      stat_score[k]=result[k][l];
    }
    score=average_score(stat_score);
    if (l<NFUNCTIONS-1)
      printf("%f\n",score);
    else
       printf("%f",score);
  }
}
