/* Function to compute diffusion of ligands.Geometry encodes the neighbouring cells, i.e. 
geometry[i][] contains the indexes of the neighbours of cell i.  The table of diffusion
csts, difflig[] is static global variable in header
*/


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
double N_to_TAU(int);
int N_to_L(int);
double N_to_L_SELF(int);




 
 void print_history( int try )  {
 
     int pas, i, j, k;
     char titre[50];
     FILE *fileptr;
 
     sprintf(titre, "Buffer%i", try);  	//the number appended to the file Buffer is the "trial" number.
     fileptr=fopen(titre, "w");
 //     for(k=0;k<NOUTPUT;k++){
//          fprintf(fileptr,"%d\n",trackout[k]);
//      }
//     fprintf(fileptr,"11235813213455\n");
fprintf(fileptr,"%i", 1);
			for(j=0; j<NCELLTOT; j++)  {
			    //fprintf(fileptr,"%2.0f\t%5d\t%.0f",N_to_TAU(j),N_to_L(j),N_to_L_SELF(j));
	    		for(i=0; i<SIZE; i++)  {
		    		fprintf(fileptr,"\t%10.15f",history[i][1][j]);
	    		}
	    		//fprintf(fileptr,"\n");
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





double compute_noise_increment(double rate, double dt)
{/*computes the increment to add to a ds*/

return gaussdev()*sqrt(rate/(dt*CONCENTRATION_SCALE));


}
