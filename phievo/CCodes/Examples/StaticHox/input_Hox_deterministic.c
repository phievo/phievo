/* define input variables as a function of time.  n-attempts = loop [0,ntries)  */



static double variability=0;


void inputs(int pas,int ncell, int ntry){
int trackg1=trackin[0];
if ((pas==0)&&(ncell==0))
   {      variability=0.1*FRAND();
   }



 double ss=exp(0.07*ncell)/8.0*(1+variability);

  
 history[trackin[0]][pas][ncell]=ss;
 if (pas>NSTEP/2)
   history[trackin[0]][pas][ncell]=ss*exp(-0.1*(pas-NSTEP/2)*DT);
    
}


