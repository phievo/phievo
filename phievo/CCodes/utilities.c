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






/* print history array to file= BUFFER, where int trial is 0,1,..NTRIES-1 */

void print_history( int trial )  {

    int pas, i, j;
    char titre[50];
    FILE *fileptr;

    sprintf(titre, "Buffer%i", trial);
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
