/* define input variables as a function of time.  n-attempts = loop [0,ntries) */

static int init_index=0;
void inputs(int pas,int ncell,int n_attempts){
  //if ((pas==0)&&(ncell==0))
  //  PHI=FRAND()*6.2832;
  //if ((pas==NSTEP/6)&&(ncell==0))//random phase shift
    // PHI=FRAND()*6.2832;


    static int ncalls;
 
  int trackg1=trackin[0];
 
    double RT=pas*DT;//Real time
   
    
      
    double dummy= 1+cos(6.28*RT/TIME*n_cycle+PHI);
    int l;
    for (l=1;l<10;l++){
      if ((pas>l*NSTEP/10)&&(ncell>0))
	dummy=(ncell*0.1*l*0.1+0.2); //Input constant from cycle 4 for the other cases 
      /*if ((pas>2*NSTEP/3+ncell*15)&&(ncell>0))
	dummy=1.5;*/
    }
  
    
    
   history[trackg1][pas][ncell]=dummy;
	  
 
}


