/* define input variables as a function of time.  n-attempts = loop [0,ntries) */

void inputs(int pas,int ncell,int n_attempts){
  int trackg1=trackin[0];
  double ss;
 
     ss=ncell-0.6*pas*DT;
     if (ss>0) 
     	ss=1;
     else 
     	ss=0;
     
     history[trackg1][pas][ncell]=ss;

}


