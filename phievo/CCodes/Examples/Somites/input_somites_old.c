/* define input variables as a function of time.  n-attempts = loop [0,ntries) */

void inputs(int pas,int ncell,int n_attempts){
  int trackg1=trackin[0];
  double ss;
 
     ss=20*exp(0.15*ncell-0.6*pas*DT);
     if (ss>500) ss=500;
     history[trackg1][pas][ncell]=ss;

}


