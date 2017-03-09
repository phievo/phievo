 /* Create time dependent input from [0, NSTEP) by integrating a sum of random positive delta functions
and then exponentiating.  The fitness demands that output after such a pulse is at least 0.2 a.u.  higher that output before the pulse
*/

static double isignal[NSTEP][NCELLTOT];
static double dsignal[NSTEP][NCELLTOT];

#define NFUNCTIONS 3 //number of  functions computed by the fitness function. should be at least 1 for the fitness

static double result[NTRIES][NFUNCTIONS]; //global matrix containing the all the various fitnesses

    
static int t1 = 1;   // duration of delta function spike in dsignal, dt units
static int t2 = 500;  // mean time of pulse

static int tfitness=NSTEP/6;//time from which we start pulses of I and compute the fitness

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

	tnext += tfitness+(int) t2*(rand()/( (double) RAND_MAX +1));
	ampl = 5+10*(rand()/( (double) RAND_MAX +1));

	for( t=tnext; t<tnext+t1; t++) {
	  dsignal[t][k] = ampl;
	}
	tnext += t1;
	
	tnext += 100+(int) 2*t2*(rand()/( (double) RAND_MAX +1));
	for( t=tnext; t<tnext+t1; t++) {
	  dsignal[t][k] -= ampl;
	}
	
    }

	
/* integrate the derivative and expon to get real signal applied to cell  */

    for(k=0; k<NCELLTOT; k++)  {
      isignal[0][k] = 0;
      for(t=1; t<NSTEP; t++) {
	isignal[t][k] = isignal[t-1][k] + DT*dsignal[t][k];
      }
      for(t=0; t<NSTEP; t++) {
	isignal[t][k] = exp(isignal[t][k])-1;
      }
    }

}

void nogood(int ntry){
  //dummy function to update result if the fitness is no good
    result[ntry][0]=RAND_MAX;
    result[ntry][1]=RAND_MAX;
    result[ntry][2]=0.0;
 
}
    
void  fitness( double history[][NSTEP][NCELLTOT], int trackout[],int ntry)  {

    static int ncalls;
    static double sum_scores;
    static int max_lag = 100;
    double norm, *hist1;
    int ll;
    int k,t;
    int delay=10;
    int ntrans=0;
    norm =0;
    double av=0;
    int i=0;
    for(k=0; k<NCELLTOT; k++)  {
      for (i=0;i<SIZE;i++){
	if (history[i][NSTEP-1][k]<0)
	  {
	    nogood(ntry); return;
	    }
	if (history[i][NSTEP-1][k]<1e10)
	  {}
	else
	  {nogood(ntry); return;
	   //intercepts nan
	  }
      }
    }
    int tinduction=0;
    int tfall=0;
    int indexinduction=0;
    int indexfall=0;
    int tprevious=tfitness;
    double diff1, diff2,diff3,currdiff1,currdiff2;
    double olddiff1, olddiff2;
   
    int ngene=0;
    for(k=0; k<NCELLTOT; k++)  {
      
      olddiff1 = olddiff2 = 0;

      if (fabs(history[trackout[0]][tfitness/2][k]/history[trackout[0]][tfitness][k]-1)>0.05)
	{
	  nogood(ntry); return;
	  //if not within 5% of steady state before the first pulse, does not select
	}
      
      

      
      for(t=tfitness; t<NSTEP; t++){
	if (history[trackin[0]][t][k]!=history[trackin[0]][t-1][k]){ //detects the pulse
	  if (tinduction==0)
	    tinduction=t;
	  else
	    tfall=t; 
	}
      }

      //make sure steady state is reached at the end
      if( fabs(history[trackout[0]][NSTEP-1][k]/history[trackout[0]][(NSTEP-1+tfall)/2][k]-1) > 0.05){
	nogood(ntry); return;
      }	

      diff1 = fabs(history[trackout[0]][tfall-1][k]-history[trackout[0]][tfitness][k]);
      diff2 = fabs(history[trackout[0]][tfall-1][k]-history[trackout[0]][NSTEP-1][k]);
      diff3 = fabs(history[trackout[0]][tfitness][k]-history[trackout[0]][NSTEP-1][k]);
      // absolute change must be at least 0.1 and relative change must be at least 20%
     //must return to the same point
      if( diff1 < 0.1 || diff1/history[trackout[0]][tfitness][k] < 0.2 || diff3/history[trackout[0]][tfitness][k]){ 
	nogood(ntry); return;
      } 
	// for on accept any change, for off only approach to steady state

//for on time, compute time to change by 5%
   for(t=tinduction; t<NSTEP; t++){
	if( t < tfall && indexinduction == 0){
        	currdiff1 = fabs(history[trackout[0]][t][k]-history[trackout[0]][tfitness][k]); 
                if( currdiff1 > 0.05 ){ indexinduction = t; break;}
	}
    }

// for off time, work backwards from steady-state, find final approach to within 5% of steady state
	for( t = NSTEP-1; t > tfall; t--){
        	currdiff2 = fabs(history[trackout[0]][t][k]-history[trackout[0]][NSTEP-1][k])/history[trackout[0]][NSTEP-1][k]; 
		 if( currdiff2 > 0.05  && indexfall == 0){
	  		indexfall=t; break;
		 }
	}
  
    //averages over the pulses the values of the max deviation and of the final deviation 
    //printf(" %i %i %i %i\n",indexinduction,tinduction,indexfall,tfall);
    int ton=indexinduction-tinduction;
    int toff=indexfall-tfall;
    if (ton>0)
      result[ntry][0]=toff/(1.0*ton);
    result[ntry][1]=toff;
    result[ntry][2]=ton;
    if (ton<=0)
      {result[ntry][0]= RAND_MAX;
	result[ntry][2]=0;
      }
    if (toff<=0)
      {result[ntry][0]= RAND_MAX;
	result[ntry][1]= RAND_MAX;
      }
  
}
}
void treatment_fitness( double history2[][NSTEP][NCELLTOT], int trackout[]){
 
  int l,k;
  double score_on,score_off,score;
  double stat_score[NTRIES];
  double stat_on[NTRIES];
  double stat_off[NTRIES];
  double average_on, average_off;
  int count_on, count_off, on_bad, off_bad;

  average_on = average_off = 0.0;
  count_on = count_off = on_bad= off_bad = 0;
  for (k=0;k<NTRIES;k++){
    stat_off[k]=1.0*result[k][1];
    stat_on[k]=1.0*result[k][2];
    stat_score[k]=1.0*result[k][0];
  }

for( k = 0 ; k < NTRIES; k++){
	if( stat_on[k] < 10000){
		average_on+=stat_on[k];
		count_on++;
	}
	else{
		on_bad++;
	}

	if( stat_off[k] < 10000 ){
		average_off+=stat_off[k];
		count_off++;
	}
	else{
		off_bad++;
	}
	
}
	if( count_on > 0)
		average_on=average_on/count_on;
	else
		average_on = 0;

	if( count_off > 0)
		average_off=average_off/count_off;
	else
		average_off=10000;

	score_on=average_on*(1-1.0/NTRIES*on_bad);
	score_off=average_off*(1+1.0/NTRIES*off_bad);
	score=score_on/score_off;
 


      printf("%f\n%f\n%f",-score_on,score_off,score);

      

  /*for (l=1;l<NFUNCTIONS;l++){
    for (k=0; k<NTRIES; k++){
      stat_score[k]=result[k][l];

    }
    score=average_score(stat_score);
    if (l<NFUNCTIONS-1)
      printf("%f\n",score);
    else
       printf("%f",score);
       }*/
}
