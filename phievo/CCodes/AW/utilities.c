/* Function to compute diffusion of ligands.Geometry encodes the neighbouring cells, i.e. 
geometry[i][] contains the indexes of the neighbours of cell i.  The table of diffusion
csts, difflig[] is static global variable in header
*/

void diffusion(int ncell, int step, double ds[],double history[][NSTEP][NCELLTOT], int geometry[][NNEIGHBOR]){

  int g,neig,index,index_ligand,n_localneig;
  double diff,diffusion;
 
 
  for (g=0;g<NLIGAND;g++){
    if (externallig[g]==1){//computes diffusion only for external ligands
      diffusion=0;
      index_ligand=tracklig[g];//tracks the ligand
      diff=difflig[g];//takes the corresponding diffusion constant
      n_localneig=0;//computes number of local neighbours
      for (neig=0;neig<NNEIGHBOR;neig++){
	index=geometry[ncell][neig];//takes the neighoubring cell
	if (index>=0){
	  diffusion+=history[index_ligand][step][index];//concentration of the ligand in the neighbouring cell
	  n_localneig+=1;
	}
      }
      diffusion-=n_localneig*history[index_ligand][step][ncell];//minus number of local neighbours times concentration of the ligand in the local cell
      diffusion*=diff;// times diffusion constant
      ds[index_ligand]+=diffusion;
    }
  }

}

/* 
computation of the ligand concentration seen by each cells
*/


/*



void sum_concentration(int ncell, int step, double concentrations[]) {

  int g,neig,index;
  // we add the ligands concentrations of all neighbouring cells 
  for (g=0;g<SIZE;g++){
    concentrations[g]=0;
    for (neig=0;neig<NNEIGHBOR;neig++){
      index=geometry[ncell][neig];
      if ((index>=0) && (index!=ncell) ) {
	concentrations[g]+=history[g][step][index];
      }
    }
  }
  // For external ligands, we simply take the local external value as local ligand concentration

  for (g=0;g<NLIGAND;g++){
    if (externallig[g]==1){
      concentrations[tracklig[g]]=history[tracklig[g]][step][ncell];
    }

  }
  



}

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

