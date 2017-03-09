/* define input variables as a function of time.  n-attempts = loop [0,ntries)  */



static double variability=0;


void inputs(int pas,int ncell, int ntry){
int trackg1=trackin[0];
if ((pas==0)&&(ncell==0))
   {variability=0.2*FRAND();
   }

 variability=0;
 double ss=exp(1.5*(ncell+1)-0.75*(1+variability)*(pas-10)*DT);
 if (pas<10)
  ss=1;
 if (ss>=1)
   ss=1;
   else
     ss=0;
  
 history[trackin[0]][pas][ncell]=ss;
 
    
}


