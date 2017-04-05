/* define input variables as a function of time.  n-attempts = loop [0,ntries) */

void inputs(int pas,int ncell,int n_attempts){
  int trackg1=trackin[0];
  double ss;
 
     ss=ncell*0.25/NCELLTOT;
     
     history[trackg1][pas][ncell]=ss;
     
     if ((pas==0)&&(n_attempts==1))
     	history[trackout[0]][pas][ncell]=10;

}


