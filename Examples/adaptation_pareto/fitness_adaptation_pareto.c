*

  Defines fitness function
*/

#define NFUNCTIONS 5 //number of  functions computed by the fitness function. should be at least 1 for the fitness

static double result[NTRIES][NFUNCTIONS]; //result will contain the fitness plus the other results computed by the fitness function

static int tfitness=NSTEP/6;//time from which we start pulses of I and compute the fitness



void nogood(int ntry){
  //dummy function to return infinite results if some fitness criteria - like crazy concentrations- is not realized
  result[ntry][0]=RAND_MAX;
  result[ntry][1]=RAND_MAX;
  result[ntry][2]=0.0;
  result[ntry][3]=RAND_MAX;
  result[ntry][4]=RAND_MAX;
}
    

void fitness( double history[][NSTEP][NCELLTOT], int trackout[],int ntry)  {

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
	}
      if (history[i][NSTEP-1][k]<1e10)
	{}
      else
	{nogood(ntry);
	  ;//intercepts nan
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
      
    int tlimit=(int) (tfitness/2);
    if (fabs(history[trackout[0]][tfitness/2][k]/history[trackout[0]][tfitness][k]-1)>0.05)
      {
	nogood(ntry);
	//if not within 5% of steady state before the first pulse, does not select
	return;
      }
      
      
    if (fabs(history[trackout[0]][tfitness][k])<0.2)
      {
	nogood(ntry);
	//minimum concentration
	return;
      }
      

      
    for(t=tfitness; t<NSTEP; t++){
       
      for(ngene=0;ngene<SIZE;ngene++){
	som_tot+=history[ngene][t][k];
      }

      som_tot-=history[trackin[0]][t][k];
	

      if (history[trackin[0]][t][k]!=history[trackin[0]][t-1][k]){ //detects the pulse of input
	  
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
      Oaverage+=history[trackout[0]][t][k];//updates average output
      Omax=MAX(Omax,history[trackout[0]][t][k]);
      Omin=MIN(Omin,history[trackout[0]][t][k]);	
    }    
  }

  
  //averages over the pulses the values of the max deviation and of the final deviation 
  sum_deviation=sum_deviation/(notfirst-1);
  sum_deriv=sum_deriv/(notfirst-1);
  result[ntry][0]=(sum_deviation+0.01)/(sum_deriv+0.0001);
  result[ntry][1]=sum_deviation;
  result[ntry][2]=sum_deriv;
  result[ntry][3]=som_tot/(NSTEP-tfitness);//average number of proteins
  result[ntry][4]=Oaverage/(NSTEP-tfitness);//average number of output
  //computes fitness
  if ((result[ntry][1]>1000)||(result[ntry][2]>1000))
    nogood(ntry);

}


void treatment_fitness( double history2[][NSTEP][NCELLTOT], int trackout[]){
 
  /* function to print out the result*/

  //if you want to do anything with the average output history2, this is the right place


  int l,k;
  double score;
  double stat_score[NTRIES];
 

  for (l=0;l<NFUNCTIONS;l++){
    for (k=0; k<NTRIES; k++){
      stat_score[k]=result[k][l];
    }
    //computes and prints the average values of the fitnesses
    score=average_score(stat_score);
    if (l<NFUNCTIONS-1)
      printf("%f\n",score);
    else
      printf("%f",score);
  }
}

