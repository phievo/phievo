 /* General information based fitness score (always returns a negative number, more neg more fit), correspond to old one with mutual information between position and states
  

*/




static double result[NTRIES];


static int t=10;


static   double TIME=NSTEP*DT;//Total time
static int n_cycle=24*NSTEP*DT/(3000*0.05);
static  double PHI= 0.0;
static  double AMP=1;//amplitude of the cosine function
 


void fitness( double history[][NSTEP][NCELLTOT], int trackout[],int ntry)  {
   
  double fitness=0;
  result[ntry]=0;
    double I2=0;
    double O2=0;
    double IO=0;
    double Iav=0;
    double Oav=0;
    double correl=0;
    int t,k;
 



for (k=0;k<NCELLTOT;k++)
  {


      for(t=0; t<NSTEP; t++){	
	double I=0;
	if (k>0)
	  I=history[trackout[0]][t][0];//reference signal in cell 0
	else
	  {//for k==0, entrainment, reference signal is the input

	    //if (cos(6.28*t*DT/TIME*n_cycle+PHI)>0)
	      I=history[trackin[0]][t][0];

	  }

	double O=history[trackout[0]][t][k];//compensated signal
	I2+=I*I;
	O2+=O*O;
	IO+=I*O;
	Iav+=I;
	Oav+=O;
      }


 I2/=NSTEP;
 O2/=NSTEP;
 IO/=NSTEP;
 Iav/=NSTEP;
 Oav/=NSTEP;
 correl=fabs(IO-Iav*Oav)/(sqrt(O2-Oav*Oav)*sqrt(I2-Iav*Iav));
 result[ntry]+=correl; //impose that the compensated signals are the same as the entrained signal

}


  result[ntry]=1-result[ntry]/NCELLTOT;
  
  
  
}

void treatment_fitness( double history2[][NSTEP][NCELLTOT], int trackout[]){
 
  /* function to print out the result*/
  //fitness(history2,trackout,NTRIES);
  double average=0;
  int k;
  for (k=0;k<NTRIES;k++)
    average+=result[k]/NTRIES;

  

  //if you want to do anything with the average output history2, this is the right place




  double zero = 0.;
  if (average<100)
    printf("%f\n%f\n%f", average, zero, zero );
  else
     printf("100\n%f\n%f", zero, zero );
  return;

}

