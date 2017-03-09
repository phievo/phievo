 /* Create time dependent input from [0, NSTEP) by integrating a sum of random positive delta functions
and then exponentiating.  The fitness demands that output after such a pulse is at least 0.2 a.u.  higher that output before the pulse
*/

static double isignal[NSTEP][NCELLTOT];
static double dsignal[NSTEP][NCELLTOT];

#define NFUNCTIONS 3 //number of  functions computed by the fitness function. should be at least 1 for the fitness

static double result[NTRIES][NFUNCTIONS]; //global matrix containing the all the various fitnesses

    
static int t1 = 1;   // duration of delta function spike in dsignal, dt units
static int t2 = 500;  // mean time of pulse

static int tfitness=NSTEP/2;//time of pulse

void init_signal( ) {

    int k, t, tnext;
    double ampl;	// random value of dsignal

    for(k=0; k<NCELLTOT; k++)  {
      for(t=0; t<NSTEP; t++) {
 		ampl = 0.1*(2.0*rand()/( (double) RAND_MAX +1)-1);
	 dsignal[t][k] = ampl;   //if multiple calls need zero expl
      }
    }

    
    tnext=tfitness;
    ampl = 5+10*(rand()/( (double) RAND_MAX +1));
    for( t=tnext; t<tnext+t1; t++) {
		dsignal[t][k] = ampl;
    }
	

	
/* integrate the derivative and expon to get real signal applied to cell  */

    for(k=0; k<NCELLTOT; k++)  {
      isignal[0][k] = (rand()/( (double) RAND_MAX +1));
      for(t=1; t<NSTEP; t++) {
	isignal[t][k] = isignal[t-1][k] + DT*dsignal[t][k];
	if (isignal[t][k] < 0) isignal[t][k] = 0;
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

    int k,i,t,num1,num2;
	double avg1,avg2,stddev1,stddev2;
	num1=num2=0;
	avg1=avg2=stddev1=stddev2=0;

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
   for(t=tfitness/2; t<tfitness; t++){
		num1++;
		avg1+=history[trackout[0]][t][k];
		stddev1+=history[trackout[0]][t][k]*history[trackout[0]][t][k];
   }
   for(t=3*tfitness/2; t<NSTEP; t++){
		num2++;
		avg2+=history[trackout[0]][t][k];
		stddev2+=history[trackout[0]][t][k]*history[trackout[0]][t][k];
   }
   avg1 = avg1/num1; avg2=avg2/num2; 
   stddev1=stddev1/num1; stddev2=stddev2/num2;
   stddev1=stddev1-avg1*avg1;
   stddev2=stddev2-avg2*avg2;

   if( stddev1 > 0 ) stddev1= sqrt(stddev1);
   else stddev1 = 0;
   if( stddev2 > 0 ) stddev2= sqrt(stddev2);
   else stddev2 = 0;

   result[k][0]=fabs(avg1-avg2);
   result[k][1]=stddev1+stddev2;

}

void treatment_fitness( double history2[][NSTEP][NCELLTOT], int trackout[]){
 
  int l,k;
  double avg_diff,avg_dev;
  double stat_dev[NTRIES];
  double stat_avg[NTRIES];
  for (k=0;k<NTRIES;k++){
    stat_dev[k]=result[k][1];
    stat_avg[k]=result[k][0];
  }

  
  avg_diff=average_score(stat_avg);
  avg_dev=average_score(stat_dev);
  if (avg_diff<1e5 && avg_dev < 1e5)
    {


      printf("%f\n%f\n%f",-avg_diff,avg_dev,avg_dev-avg_diff);

      
    }
    else
      printf("%i\n%i\n%i",RAND_MAX,RAND_MAX,RAND_MAX);

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
