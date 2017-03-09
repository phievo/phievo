 /* General information based fitness score (always returns a negative number, more neg more fit), correspond to old one with mutual information between position and states
  

*/

#define MAXSTATES 64

#define NFUNCTIONS 2 //number of  functions computed by the fitness function, >=1
static double result[NTRIES][NFUNCTIONS]; //global matrix with fitness fns each try

/* These two arrays of stored data used to compute prob[][]  */

static double output_ncell[NCELLTOT][NOUTPUT];  // save the vector of output variables


static int input_state=0;   // set in inputs() subroutine input_state = [0,.. NUMLEVELS)

static double min_history = 0.1;  // add rand value to history to kill info if h->0
static double random_final_step = 0.9;  // final time history = [this fract, 1)*NSTEP


void nogood(int ntry){
  //dummy function to fill result if there is a problem (nan, negative concentration...)

  int i=0;
  for (i=0;i<NFUNCTIONS;i++){
     result[ntry][i]=RAND_MAX;
  }
}

double limit_history(double h) {
  /* if invalid history value, replace with RAND_MAX  */

  if( 0 <= h && h < (double)RAND_MAX ) {
    return h;
  }
  return (double) RAND_MAX;
}




double mutual_information( int nrow, int ncol, double prob[][MAXSTATES] ) {
  /* return - sum p_ij log(p_i*p_j/p_ij) NB return a number > 0.
   Do not assume prob[][] is normalized */

  int i,j;
  double sum_prob=0, sum_rows[ncol], sum_cols[nrow], info=0;

  for( i=0; i<nrow; i++) {
    sum_cols[i] = 0;
    for( j=0; j<ncol; j++ ) {
      sum_cols[i] += prob[i][j];
    }
  }

  for( j=0; j<ncol; j++ ) {
    sum_rows[j] = 0;
    for( i=0; i<nrow; i++) {
      sum_rows[j] += prob[i][j];
    }
  }

  for( i=0; i<nrow; i++) {
    sum_prob += sum_cols[i];
  }

  for( i=0; i<nrow; i++) {
    for( j=0; j<ncol; j++ ) {
      prob[i][j] /= sum_prob;
      // printf("(in, out) %i %i %f\n",  i,j,prob[i][j] );
      if( prob[i][j] > FLT_MIN ) {
        // printf("%i %i %f %f %f %f\n", i,j,prob[i][j], sum_rows[j], sum_cols[i], sum_prob );
        info += prob[i][j]*log(sum_prob*sum_prob*prob[i][j]/(sum_rows[j]*sum_cols[i]) );
      }
    }
  }

  return info;
}








void output_prob( double prob[][MAXSTATES] ) {
  /* take output[ntries] array and compute compute un-normed prob.  
     Let p0, p1 be weights on 0 or 1 for one output concentration
     prob[][i] = p1 for output i = 0,1,2..  
     Sum over all ntry and each output by unit vector, eg out = (3,1), 
     give weight 3/4 to O1 and weight 1/4 to O2
     if no output is expressed, prob[0][NOUTPUT] is negative ->dead state
  */

  double sum, prob1, outs[NOUTPUT];
  int n,i;
 
  for(n=0; n<NCELLTOT; n++ ) {
     int track_zero=0;
    sum = 0; 

    for(i=0; i<NOUTPUT; i++) {
      outs[i] = output_ncell[n][i];
      if (outs[i]>2*min_history)
	sum += outs[i];
      else 
	outs[i]=0;
    }
   
    
   
    for(i=0; i<NOUTPUT; i++) {
      outs[i] = outs[i]/sum;
    }

    for(i=0; i<NOUTPUT; i++) {
    
	 prob[n][i] += outs[i];
	 if (sum>0.1)//if any output is high enough 
	   {	 
	     track_zero=1;
	   }
      
    }

    if (track_zero==0)//if no output is expressed, dead state
      {
	for(i=0; i<NOUTPUT; i++) {
	  prob[n][i]=0;
	}

	  prob[0][NOUTPUT]=-1.0;
	
       
      }


  }
} 








void compute_entropies(double *fitness1, double *fitness2, int nrow, int ncol, double prob[][MAXSTATES] ) {
    /* return - sum p_ij log(p_i*p_j/p_ij) NB return a number > 0.
     Do not assume prob[][] is normalized */
    
    int i,j;
    double sum_prob=0, sum_rows[ncol], sum_cols[nrow], info=0, pj, entropy=0;
    
    for( i=0; i<nrow; i++) {
        sum_cols[i] = 0;
        for( j=0; j<ncol; j++ ) {
            sum_cols[i] += prob[i][j];
        }
    }
    
    for( j=0; j<ncol; j++ ) {
        sum_rows[j] = 0;
        for( i=0; i<nrow; i++) {
            sum_rows[j] += prob[i][j];
        }
    }
    
    for( i=0; i<nrow; i++) {
        sum_prob += sum_cols[i];
    }
    
    for( i=0; i<nrow; i++) {
        for( j=0; j<ncol; j++ ) {
            prob[i][j] /= sum_prob;
            // printf("(in, out) %i %i %f\n",  i,j,prob[i][j] );
            if( prob[i][j] > FLT_MIN ) {
                // printf("%i %i %f %f %f %f\n", i,j,prob[i][j], sum_rows[j], sum_cols[i], sum_prob );
                info += prob[i][j]*log(sum_prob*sum_prob*prob[i][j]/(sum_rows[j]*sum_cols[i]) );
            }
           
        }

    }
    
    for( j=0; j<ncol; j++ ) {
        pj=sum_rows[j]/sum_prob;
        if (pj> FLT_MIN)
        {
            entropy-=pj*log(pj);
        }
    }
    
    if (prob[0][NOUTPUT]>=0){
        *fitness1=entropy;
        *fitness2=entropy-info;
    }
    else
    {
        *fitness1=-1.0;
        *fitness2=100.0;
    }
}







void output2fitness( double *fitness1, double *fitness2) {
  /* From the vector of output for each try, compute a probability, then implements entropies/mutual info  */

    double prob[NCELLTOT][MAXSTATES] = {0};
   
    output_prob( prob );
    compute_entropies(fitness1,fitness2, NCELLTOT,NOUTPUT,prob);

}

void fitness( double history[][NSTEP][NCELLTOT], int trackout[],int ntry)  {
     int ncell, ngene, nout, last_step;
  double random_final_state = 0.8;
  double fitness1, fitness2, final_concentration[SIZE][NCELLTOT] = {0}, tmp;

  /* Computes the fitness of individual embryos*/
  // randomize the last step and add up all cells
  tmp = ( random_final_state + (1 - random_final_state)*FRAND() )*NSTEP;
  last_step = (int) tmp;
  for (ncell=0; ncell<NCELLTOT; ncell++){

    for (ngene=0; ngene<SIZE; ngene++){
    
      final_concentration[ngene][ncell] = limit_history( history[ngene][last_step][ncell] ) + min_history;

  
      
    }
  }

for (ncell=0; ncell<NCELLTOT; ncell++){
  for (ngene=0;ngene<NOUTPUT;ngene++){

    output_ncell[ncell][ngene]=final_concentration[trackout[ngene]][ncell];
  }
 }
  
 
  
  for (ncell=0; ncell<NCELLTOT; ncell++){

    for (ngene=0; ngene<SIZE; ngene++){
     if( final_concentration[ngene][ncell] >= (double) RAND_MAX ) {
    printf("%f\n",(double) RAND_MAX );
  }
  
      
    }
  }

 
  
  output2fitness( &fitness1, &fitness2);




  result[ntry][0]=fitness1;
  result[ntry][1]=fitness2;

  
}

void treatment_fitness( double history2[][NSTEP][NCELLTOT], int trackout[]){
 
  /* function to print out the result*/



  //if you want to do anything with the average output history2, this is the right place

  /* For this specific problem, we actually compute the fitness only on the "averaged" embryo" , we do not use the table result*/
    
  int ncell, ngene, nout, last_step,ntry;
  double random_final_state = 0.9;
  double fitness1, fitness2,final_concentration[SIZE][NCELLTOT] = {0}, tmp;

  
  // randomize the last step and add up all cells
  tmp = ( random_final_state + (1 - random_final_state)*FRAND() )*NSTEP;
  last_step = (int) tmp;
  for (ncell=0; ncell<NCELLTOT; ncell++){

    for (ngene=0; ngene<SIZE; ngene++){
    
      final_concentration[ngene][ncell] = limit_history( history2[ngene][last_step][ncell] ) + min_history;

  
      
    }
  }

for (ncell=0; ncell<NCELLTOT; ncell++){
  for (ngene=0;ngene<NOUTPUT;ngene++){

    output_ncell[ncell][ngene]=final_concentration[trackout[ngene]][ncell];
  }
 }
  
 
  
  for (ncell=0; ncell<NCELLTOT; ncell++){

    for (ngene=0; ngene<SIZE; ngene++){
     if( final_concentration[ngene][ncell] >= (double) RAND_MAX ) {
    printf("%f\n",(double) RAND_MAX );
  }
  
      
    }
  }

 
  
  output2fitness( &fitness1, &fitness2);


   
  




 /* penality if no steady state   */
    double penal=0;
  for (ngene=0;ngene<NOUTPUT;ngene++)  {
    for (ncell=0;ncell<NCELLTOT;ncell++)  {
        int n_pas;
        for (n_pas=NSTEP-100;n_pas<NSTEP;n_pas+=5)  {
	  tmp = limit_history( history2[trackout[ngene]][NSTEP-1][ncell] );
	  if( tmp == (double) RAND_MAX ) {
	    penal += tmp;
	  }
	  else {
	    penal += fabs(limit_history( history2[trackout[ngene]][n_pas][ncell] )- tmp );
	  }
        }
     }
  }



  double zero = 0.;
  if (penal<100) //works with 0.1, trying 0.01
    printf("%f\n%f", -fitness1+0.01*FRAND(), fitness2+0.01*FRAND() );
  else
     printf("%f\n%f", zero, 100.0 );
  return;




}

