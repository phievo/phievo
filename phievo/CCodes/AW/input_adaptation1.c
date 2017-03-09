/* define input variables as a function of step, which cell, and n-attempts = loop [0,ntries)
set history[input# ][pas][ncell] to correct value.
*/

void inputs(int pas,int ncell,int n_attempts){

   int track = trackin[0];

   history[track][pas][ncell] = isignal[pas][ncell];
}


