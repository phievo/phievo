/* Template for the fitness, the fitness() and treatment_fitness() are used by the
main.c and must be supplied with the arguments and returns indicated.
   The fitness() called for each 'try' (eg iteration over initial conditions and inputs)
Any results needed for computing fitness for selection and other diagonstic to be
returned to python code are accumulated in result[][] array.   
   treatment_fitness is called by main after all ntries are run and 
computes what is to be returned to python.  At this stage can compute the fitness
from history[][][] arveraged over all tries.
*/

#define NFUNCTIONS 2 //number of  functions computed by the fitness function. should be at least 1 for the fitness

static double result[NTRIES][NFUNCTIONS]; //global matrix containing the all the various fitnesses

    

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

    double amax = 0;
    int ind0 = trackout[0];
   
    //find max of output 
    for( i = 0 ; i < NSTEP; i++){
	if( history[ind0][i][0] > amax ) amax = history[ind0][i][0];
    }
    // select for nonmonotonic (number between -1 and 1)
    result[ntry][0] = (history[ind0][0][0]+history[ind0][NSTEP-1][0]-amax)/
			(history[ind0][0][0]+history[ind0][NSTEP-1][0]+amax);
    
    
    double mut = 0;
    
    for( i = 0 ; i < NSTEP; i++){	
	//integrate over enzyme concentrations with saturation
	mut = mut+history[ind0][i][0]/(1+history[ind0][i][0])*DT;
    } 

    result[ntry][1] = mut; //store total mutations
    
    
  
  
}


void treatment_fitness( double history2[][NSTEP][NCELLTOT], int trackout[]){
 
  /* function to print out the result*/



  //if you want to do anything with the average output history2, this is the right place


  int l,k;
  double score,nonmon,corr,avgI,avgM,devI,devM;
  double stat_score[NTRIES];
  for (k=0;k<NTRIES;k++){
    stat_score[k]=result[k][1];
  }

  /* define score that will be used for selection here.  For most problems 
  score and average are the same  */

  nonmon = 0; 
  for( k = 0; k < NTRIES; k++){
	nonmon+=result[k][0]/NTRIES;
  }

  avgI = 5.0; // update depending on how input distribution is chosen
  devI = 3.2;
  avgM = average_score(stat_score);
  devM = std_dev(stat_score);

  corr = 0.0;
  for( k = 0 ; k < NTRIES; k++){
	corr+=(result[k][1]-avgM)*(10.0*k/9.0-avgI);
  }

  corr = corr/(NTRIES*devM*devI); //pearson correlation
  
  score = nonmon+corr;

  /*score=average_score(stat_score);
  double average=average_score(stat_score);
  double variance_score=std_dev(stat_score); */

  if (score < RAND_MAX){
      //printf("%f\n%f\n%f\n",score,average,variance_score);     
      printf("%f\n%f\n%f\n",score,nonmon,corr);     
  }
  else
      printf("%i\n%i\n%i\n",RAND_MAX,RAND_MAX,RAND_MAX);

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
