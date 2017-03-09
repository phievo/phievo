#define SIZE 2
#define NSTEP 30000
#define NCELLTOT 1
#define NNEIGHBOR 3
#define NTRIES 10
#define  CONCENTRATION_SCALE 0.000000
#define NFREE_PRMT 0
static double free_prmt[] = {}; 
#define	NOUTPUT 1
#define	NINPUT 1
#define	NLIGAND 0
#define	NDIFFUSIBLE 1
#define SEED 60206
#define PRINT_BUF 0
#define DT 0.050000
static int trackin[] = {0};
static int trackout[] = {1};
static int tracklig[] = {};
static int trackdiff[] = {1};
static double diff_constant[] = {0.567914843964};
static int externallig[] = {};


 
/********  header file, arrays and functions common to all following subroutines and
  the same for all problems.  NOTE this file preceeded by compiler directives defining
  parameters and the track*[] arrays
********/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <limits.h>

/* global arrays for history and geometry */

static double history[SIZE][NSTEP][NCELLTOT];
static int geometry[NCELLTOT][NNEIGHBOR];

double MAX(double a,double b){
	 if (a>b) return a;
 	 return b;
}
 
double FRAND()  {
	return (double) rand()/((double)RAND_MAX + 1);
}
 
double MIN( double a, double b ) {
	if (a<b) return a;
	return b;
}



double POW(double x,double n){
  return exp(n*log(x));
}





double HillR(double x,double thresh,double n)
{
	double r=exp(n*log(x/thresh));
	return 1.0/(1+r);
}
 
 
double HillA(double x,double thresh,double n)
{
	double r=exp(n*log(x/thresh));
	return r/(1+r);
 }

 /* Function to compute diffusion of ligands.Geometry encodes the neighbouring cells, i.e. 
geometry[i][] contains the indexes of the neighbours of cell i.  The table of diffusion
csts, difflig[] is static global variable in header
*/

void diffusion(int ncell, int step, double ds[],double history[][NSTEP][NCELLTOT], int geometry[][NNEIGHBOR]){

  int g,neig,index,index_diff,n_localneig;
  double diff,diffusion;
 
 
  for (g=0;g<NDIFFUSIBLE;g++){
 
      diffusion=0;
      index_diff=trackdiff[g];//tracks the diffusible species
      diff=diff_constant[g];//takes the corresponding diffusion constant
      n_localneig=0;//computes number of local neighbours
      for (neig=0;neig<NNEIGHBOR;neig++){
	index=geometry[ncell][neig];//takes the neighoubring cell
	if (index>=0){
	  diffusion+=history[index_diff][step][index];//concentration of the ligand in the neighbouring cell
	  n_localneig+=1;
	}
      }
      diffusion-=n_localneig*history[index_diff][step][ncell];//minus number of local neighbours times concentration of the ligand in the local cell
      diffusion*=diff;// times diffusion constant
      ds[index_diff]+=diffusion;
    }
  }


/* 
computation of the ligand concentration seen by each cells
*/




void sum_concentration(int ncell, int step, double concentrations[]) {

  int l,g,neig,index;
  // compute local ligands concentrations, returns the table "concentrations" 
  for (l=0;l<NLIGAND;l++){
    
  if (externallig[l]==1){
      concentrations[tracklig[l]]=history[tracklig[l]][step][ncell];   // For external ligands, we simply take the local external value as local ligand concentration
  }
  else
    {//in that case, we sum ligand concentration over all neighbouring cells
      g=tracklig[l];
      concentrations[g]=0;
      for (neig=0;neig<NNEIGHBOR;neig++){
	index=geometry[ncell][neig];
	if ((index>=0) && (index!=ncell) ) {//note that we exclude the local concentration of ligand from the sum
	  concentrations[g]+=history[g][step][index];  
	}
      }


    }

  }



}






/* print history array to file= BUFFER, where int try is 0,1,..NTRIES-1 */

void print_history( int try )  {

    int pas, i, j;
    char titre[50];
    FILE *fileptr;

    sprintf(titre, "Buffer%i", try);
    fileptr=fopen(titre, "w");

    for(pas=0; pas<NSTEP; pas++)  {
	fprintf(fileptr,"%i", pas);
	for(j=0; j<NCELLTOT; j++)  {
	    for(i=0; i<SIZE; i++)  {
		fprintf(fileptr,"\t%f",history[i][pas][j]);
	    }
	}
	fprintf(fileptr,"\n");
    }
    fprintf(fileptr,"\n");
    fclose(fileptr);
}



/* statistical tools*/

/* Function that averages the fitness scores*/
double average_score(double score[]){

  double average=0;
  int k;
  for (k=0;k<NTRIES;k++)
    average+=score[k];

  return average/NTRIES;


}


/* Function that computes the variance of fitness scores.*/
double variance_score(double score[]){

  double var=0;
  int k;
  for (k=0;k<NTRIES;k++)
    var+=score[k]*score[k];
  
  var/=NTRIES;
  double average=average_score(score);
  return var-average*average;


}



/* Function that computes the std_deviation*/
double std_dev(double score[]){

  double var=0;
  int k;
  for (k=0;k<NTRIES;k++)
    var+=score[k]*score[k];

  var/=NTRIES;
  double average=average_score(score);
  return sqrt(var-average*average);


}




double gaussdev()
{/*computes a normally distributed deviate with zero mean and unit variance
   using  Box-Muller method, Numerical Recipes C++ p293 */
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;

	if (iset == 0) {
	  do {
	    v1=2.0*FRAND()-1.0;
	    v2=2.0*FRAND()-1.0;
	    rsq=v1*v1+v2*v2;
	  } while (rsq >= 1.0 || rsq == 0.0);
	  fac=sqrt(-2.0*log(rsq)/rsq);
	  gset=v1*fac;
	  iset=1;//set flag
	  return v2*fac;//returns one and keep other for nex time
	} 
	else {
	  iset=0;
	  return gset;
	}
}





double compute_noisy_increment(double rate)
{/*computes the increment to add to a ds*/

return rate+gaussdev()*sqrt(rate/(DT*CONCENTRATION_SCALE));


}
 /***** end of header, begiining of python computed functions ***/

 void derivC(double s[],double history[][NSTEP][NCELLTOT],int step, double ds[],double memories[],int ncell){
 int index;	 for (index=0;index<SIZE;index++) ds[index]=0;//initialization
	 double increment=0;
	 double rate=0;

/**************degradation rates*****************/
	 	 rate=0.504687*s[0];
	 	 increment=rate;
	 	 ds[0]-=increment;
	 	 rate=0.310534*s[1];
	 	 increment=rate;
	 	 ds[1]-=increment;

/**************Transcription rates*****************/
 	 int k,memory=-1;
	 memory=step-0;
	 if(memory>=0){
	 	 rate=MAX(0.746213*HillA(history[0][memory][ncell],0.933635,1.845004),0.000000);
	 	 increment=rate;
	 	 ds[1]+=increment;
	}

/**************Protein protein interactions*****************/

/**************Phosphorylation*****************/
 float total;
}

void LRinC(double s[],double ds[],double ligands[]){
    double increment=0;
    double rate=0;
}
 
 /***** end of python computed functions, beginning problem specific fns ***/

  /*

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

 /* Define geometry[NCELLTOT][NNEIGHBOR] array defining the neighbors of each cell
   Convention geometry < 0 ie invalid index -> no neighbor, ie used for end cells
   NB even for NCELLTOT=1, should have NNEIGHBOR == 3 to avoid odd seg faults since
   3 neighbors accessed in various loops, eventhough does nothing.
*/
   
void init_geometry() {

  int index;

  for (index=0;index<NCELLTOT;index++){
    geometry[index][0]=index-1;//NB no left neighbour for cell 0 (-1)
    geometry[index][1]=index;
    geometry[index][2]=index+1;
  }

  geometry[NCELLTOT-1][2]=-1;//NB no right neighbour for cell ncell-1 (-1)

}


 /* set history of all species at initial time to rand number in [0.1,1) eg avoid 0 and initializes Input.
For this problem, we created a function init_signal that defines random steps of inputs.
init_signal creates time dependent input from [0, NSTEP) by integrating a sum of random positive delta functions and then exponentiating. 

tfitness is defined in the fitness C file
*/

static double isignal[NSTEP][NCELLTOT];
static double dsignal[NSTEP][NCELLTOT];

    
static int t1 = 1;   // duration of delta function spike in dsignal, dt units
static int t2 = 1000;  // mean time between delta fns (uniform distrib [1, 2*t2]  dt units


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
	tnext += 500+(int) 2*t2*(rand()/( (double) RAND_MAX +1));
	ampl = 20*(2*(rand()/((double)RAND_MAX + 1)) -1);

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


void init_history()  {

    int ncell,n_gene;
    init_signal();  // initialize signal for all times.

    for (ncell=0;ncell<NCELLTOT;ncell++){
    	for (n_gene=0;n_gene<SIZE;n_gene++){
	  history[n_gene][0][ncell] = FRAND();
        }
    }
}

 /* define input variables as a function of step, which cell, and n-attempts = loop [0,ntries)
set history[input# ][pas][ncell] to correct value.
*/


void inputs(int pas,int ncell,int n_attempts){

   int track = trackin[0];

   history[track][pas][ncell] = isignal[pas][ncell];
}


 
/* compute the RHS of equations and run over NSTEP's with 1st order Euler method
   The arugment kk, is an index passed to inputs that records which iteration of 
   initial or boundary conditions the same system of equs is being integrated for

   This version of code used with ligand at cell computed from sum neighbors
*/
 
void integrator(int kk){

    double s[SIZE];
    double ds[SIZE]; 
    double sumligands[SIZE];
    double memory[SIZE];
    int index,n,pas,ncell;

    for (index=0;index<SIZE;index++){
	s[index] = 0;
        ds[index]=0;
        memory[index]=0;
    }
 
    /* initialize geometry here, incase cells move  */
    init_geometry();
    init_history();

    /* loop over time steps, then over each cell etc */
    for (pas=0;pas<NSTEP-1;pas++)  {
	for (ncell=0;ncell<NCELLTOT;ncell++)  {
            inputs(pas,ncell,kk);
            for (index=0;index<SIZE;index++) {
	        s[index]=history[index][pas][ncell];
            }
            derivC(s,history,pas,ds,memory,ncell);  //local integration
            sum_concentration(ncell,pas,sumligands);  //perform sum of ligands concentrations for non external ligands
	    diffusion(ncell,pas,ds,history,geometry);//computes diffusion of external ligands
            LRinC(s,ds,sumligands);
	   
            for (index=0;index<SIZE;index++) {
	 	 history[index][pas+1][ncell] = s[index] + DT*ds[index];
		 if (history[index][pas+1][ncell]<0)//might happen for langevin
		   history[index][pas+1][ncell]=0;
	    }
	}
    }

    /* fill in inputs for last time.  */
    for (ncell=0;ncell<NCELLTOT;ncell++)  {
      inputs(NSTEP-1,ncell,kk);
    }
}

 

/*  General main.c, should not need any others, see comments in fitness_template for 
the functions that must be supplied in that file for the main to work,
*/

int main()  {
    srand( SEED );
    int i, k,l;
    double score = 0;
    int len_hist = SIZE*NCELLTOT*NSTEP;
    double *hptr = &history[0][0][0];
    /* following incase one wants to average history before doing fitness */
    static double history2[SIZE][NSTEP][NCELLTOT];
    double *h2ptr = &history2[0][0][0];//table for averaging output (used for multicell problems)

   

   
    /* dummy return when no outputs for fitness function */

    if(NOUTPUT <= 0) {
    	printf("%s","no output variables? terminating without integration" );
	//return;
    }
    
    for(i=0; i<len_hist; i++)  {
	*(h2ptr + i) = 0;
    }
    

    for (k=0; k<NTRIES; k++){
    	integrator(k);
	fitness(history, trackout,k);
       	if( PRINT_BUF )  {
	    print_history(k);
	}

	for(i=0; i<len_hist; i++)  {
	    *(h2ptr+i) += *(hptr +i);
	}
    }


    treatment_fitness(history2,trackout);

  
}  

