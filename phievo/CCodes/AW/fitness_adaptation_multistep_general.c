 /* Create time dependent input from [0, NSTEP) by integrating a sum of random positive delta functions
and then exponentiating.  The fitness demands that output after such a pulse is at least 0.2 a.u.  higher that output before the pulse
*/

static double isignal[NSTEP][NCELLTOT];
static double dsignal[NSTEP][NCELLTOT];

#define NFUNCTIONS 5 //number of  functions computed by the fitness function. should be at least 1 for the fitness

static double result[NTRIES][NFUNCTIONS]; //result will contain the fitness plus the other results computed by the fitness function

    
static int t1 = 1;   // duration of delta function spike in dsignal, dt units
static int t2 = 500;  // mean time between delta fns (uniform distrib [1, 2*t2+1]  dt units

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

      while( tnext < NSTEP ) {
	tnext += 1+(int) 2*t2*(rand()/( (double) RAND_MAX +1));
	ampl = 25*(2*(rand()/((double)RAND_MAX + 1))-1);

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

void nogood(int ntry){
  //dummy function to update result if the fitness is no good
    result[ntry][0]=RAND_MAX;
    result[ntry][1]=RAND_MAX;
    result[ntry][2]=0.0;
    result[ntry][3]=RAND_MAX;
    result[ntry][4]=RAND_MAX;
}
    

void fitness( double history[][NSTEP][NCELLTOT], int trackout[], int ntry)  {

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
	    nogood(ntry);
	    return;
	  }
	if (history[i][NSTEP-1][k]<1e10)
	  {}
	else
	  {nogood(ntry);
	   return;
	  }
      }
    }
    double sum_deviation=0;
    double sum_deriv=0;
    double Omin=RAND_MAX;
    double Omax=0;
    int tprevious=tfitness;
    double Oaverage=0;
    double DeltaO=0;
    double Oprevious=0;
    int notfirst=0;//we want to compute the fitness only after the second pulse
    sum_scores=0;
    int ngene=0;
    double som_tot=0;
    for(k=0; k<NCELLTOT; k++)  {
      

	/*
      int tlimit=(int) (FRAND()*(tfitness)*0.25 +tfitness/2);
      if (fabs(history[trackout[0]][tfitness/2][k]/history[trackout[0]][tfitness][k]-1)>0.05)
	{
	  nogood(ntry);
	  return RAND_MAX;//if not within 5% of steady state before the first pulse, does not select
	}
      */
      
        int tlimit=0;
	for (tlimit=tfitness/2;tlimit<tfitness;tlimit++){ 
      if (fabs(history[trackout[0]][tlimit][k]/history[trackout[0]][tfitness][k]-1)>0.05)
	{
	  nogood(ntry); return;
	  //if not within 5% of steady state before the first pulse, does not select
	}
	  
	}


      
      for(t=tfitness; t<NSTEP; t++){
       
	for(ngene=0;ngene<SIZE;ngene++){
	  som_tot+=history[ngene][t][k];
	}

	som_tot-=history[trackin[0]][t][k];
	

	if (history[trackin[0]][t][k]!=history[trackin[0]][t-1][k]){ //detects the pulse
	  
	  DeltaO=MAX(fabs(Omax-Oprevious),fabs(Omin-Oprevious));//computes maximum deviation from previous concentration of O
	  
	  if (notfirst>0)
	    {sum_deviation+=fabs(history[trackout[0]][t-1][k]-Oprevious);//computes final deviation with previous value of O 
	     sum_deriv+=DeltaO;
	   
	    }
	  
	  notfirst+=1;
	  Oprevious=history[trackout[0]][t-1][k];//updates previous value of O	  
	  //reinitializes other variables
	  Omin=RAND_MAX;
	  Omax=0;
	  tprevious=t;
	}
	Oaverage+=history[trackout[0]][t][k];
	Omax=MAX(Omax,history[trackout[0]][t][k]);
	Omin=MIN(Omin,history[trackout[0]][t][k]);


	

      }
      
    }

  
    //averages over the pulses the values of the max deviation and of the final deviation 
    sum_deviation=sum_deviation/(notfirst-1);
    sum_deriv=sum_deriv/(notfirst-1);
    
    result[ntry][0]=sum_deviation;
    result[ntry][1]=-sum_deriv/(1+sum_deriv);
    result[ntry][2]=som_tot/(NSTEP-tfitness);//average number of proteins
    result[ntry][3]=Oaverage/(NSTEP-tfitness);//average number of output
    result[ntry][4]=result[ntry][0]+result[ntry][1];

}

void treatment_fitness( double history2[][NSTEP][NCELLTOT], int trackout[]){
 
  /* function to print out the result*/



  //if you want to do anything with the average output history2, this is the right place


  int l,k;
  double score,sum_deviation,sum_deriv;
  double stat_score1[NTRIES], stat_score2[NTRIES];
  for (k=0;k<NTRIES;k++){
    stat_score1[k]=result[k][0];
    stat_score2[k]=result[k][1];
  }

  /* define score that will be used for selection here.  For most problems 
  score and average are the same  */
  
  sum_deviation = average_score(stat_score1);
  sum_deriv = average_score(stat_score2);
  score = sum_deriv/sum_deviation;
  
  

  /*score=average_score(stat_score);
  double average=average_score(stat_score);
  double variance_score=std_dev(stat_score); */

  //if (score < RAND_MAX){
      //printf("%f\n%f\n%f\n",score,average,variance_score);     
      printf("%f\n%f\n%f",sum_deviation,sum_deriv,score);  
      if( NFUNCTIONS > 1) printf("\n");   
  //}
  //else
 //     printf("%i\n%i\n%i\n",RAND_MAX,RAND_MAX,RAND_MAX);

  /*for (l=1;l<NFUNCTIONS;l++){
    for (k=0; k<NTRIES; k++){
      stat_score[k]=result[k][l];
    }
    score=average_score(stat_score);
    if (l<NFUNCTIONS-1)
      printf("%f\n",score);
    else
       printf("%f",score);
  } */
}


