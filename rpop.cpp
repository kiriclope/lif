#include "librairies.h"
#include "GlobalVars.h"
#include "Net_Utils.h"
#include "Matrix_Utils.h"

clock_t t1=clock();

double TL(double x) {  
  if(x>0.)
    return x ;
  else 
    return 0. ;
}

double F(double x) {  
  return .5 * ( 1.0 + erf( x / sqrt(2.0) ) ) ; 
}

int main(int argc , char** argv) {

  ///////////////////////////////////////////////////////////////////    
  // Parameters 
  ///////////////////////////////////////////////////////////////////    

  string dir ;
  int nbpop, ndI;
  unsigned long N ;
  double K, **J, **Tsyn, g, *Iext, *IextBL ;
  
  Set_Parameters(argc,argv,dir,nbpop,N,K,g,Iext,IextBL) ; 

  ///////////////////////////////////////////////////////////////////    

  double *Tr = new double [nbpop]() ;
  for(int i=0;i<nbpop;i++)
    Tr[i] = 10. ;

  Import_Synaptic_Parameters(nbpop,Tsyn,dir) ;

  double **Cstsyn = new double *[nbpop]() ;
  for(int i=0;i<nbpop;i++) {
    Cstsyn[i] = new double [nbpop]() ;    
    for(int j=0;j<nbpop;j++)
      Cstsyn[i][j] = exp(-dt/Tsyn[i][j]) ;
  }

  ///////////////////////////////////////////////////////////////////    
  // Scaling
  ///////////////////////////////////////////////////////////////////    

  cout << "External Inputs : " ;
  for(int i=0;i<nbpop;i++) {
    cout << Iext[i] << " " ;
    Iext[i] = sqrt(K)*Iext[i]*m0 ;
  }
  cout << endl ;

  cout << "Connectivity Parameters :" << endl ;
  Import_Connectivity_Parameters(nbpop,J,dir) ;
  for(int i=0;i<nbpop;i++) {
    for(int j=0;j<nbpop;j++) {
      cout << J[i][j] << " ";
      J[i][j] = g*J[i][j]/sqrt(K);
    }
    cout << endl ;
  }

  ///////////////////////////////////////////////////////////////////    

  double (*Func)(double) = NULL ;
  // Func = &F ;
  Func = &TL ;

  string path = "../" ;
  CreateDir(dir, nbpop, N, K, g, path) ; 

  ///////////////////////////////////////////////////////////////////    


  int i,j ;
  double tw=0,tc=0,percentage=0 ;
  double RK1=0,RK2=0,RK3=0,RK4=0,RK5=0,RK6=0;

  int IF_TRANSIENT = 0 ;
  
  ///////////////////////////////////////////////////////////////////    

  vector<double> r(nbpop) ;

  vector<double> Mean_Activity(nbpop) ; // population's rates
  vector<double> Isyntot(nbpop) ;
  vector<vector<double> > Isyn(nbpop,vector<double> (nbpop));

  ///////////////////////////////////////////////////////////////////    

  string rStr = path + "/PopRates.txt" ;
  ofstream rFile(rStr.c_str(), ios::out | ios::ate);

  ///////////////////////////////////////////////////////////////////    

  cout <<"Initialization of the network" <<endl;
  
  //Using C++11
  // default_random_engine gen  ;
  random_device rd;
  default_random_engine gen( rd() );
  uniform_real_distribution<double> unif(0,1.) ;
  normal_distribution<double> gaussian(1,.25) ;
  
  for(i=0;i<nbpop;i++) {
    r[i] = unif(gen) ; 
    for(j=0;j<nbpop;j++) 
      Isyn[i][j] = J[i][j]*r[j] ; 
  }
  
  ///////////////////////////////////////////////////////////////////    

  cout << "Main loop :" ;
  cout << " duration " << duration << " | dt " << dt ;
  cout << " | Tst " << Tst << " | Tw " << Tw << endl ;
  for (double t=0.; t<=duration; t+=dt) {
    
    percentage = t/duration ; 
    
    ///////////////////////////////////////////////////////////////////    
    // Rate dynamics
    ///////////////////////////////////////////////////////////////////    

    // for(i=0;i<nbpop;i++) {
    //   Isyntot[i] = Iext[i] ;
    //   for(j=0;j<nbpop;j++) 
    // 	if(J[i][j]!=0)
    // 	  Isyntot[i] += K * J[i][j] * r[j] ;
    // }
    
    // for(i=0;i<nbpop;i++) {
      
    //   r[i] += dt * ( -r[i] + Func(Isyntot[i]) ) / Tr[i] ;

    //   if(t>=Tst)
    // 	Mean_Activity[i] += r[i] ;
    // }
    
    ///////////////////////////////////////////////////////////////////    
    // Synaptic dynamic
    ///////////////////////////////////////////////////////////////////    

    for(i=0;i<nbpop;i++) 
      for(j=0;j<nbpop;j++) 
    	if(J[i][j]!=0)
    	  Isyn[i][j] = Cstsyn[i][j]*Isyn[i][j] + dt/Tsyn[i][j]*r[j] ;
    
    ///////////////////////////////////////////////////////////////////    
    // RK4 
    ///////////////////////////////////////////////////////////////////    
    
    for(i=0;i<nbpop;i++) { 
      for(j=0;j<nbpop;j++) 
    	if(J[i][j]!=0) { 
    	  RK1 = (-Isyn[i][j] + r[j])/Tsyn[i][j] ; 
    	  RK2 = (-Isyn[i][j]-dt/2.*RK1 + r[j] )/Tsyn[i][j] ;
    	  RK3 = (-Isyn[i][j]-dt/2.*RK2 + r[j] )/Tsyn[i][j] ;
    	  RK4 = (-Isyn[i][j]-dt*RK3 + r[j] )/Tsyn[i][j] ;
    	  Isyn[i][j] = Isyn[i][j] + dt/6.*( RK1 + 2.*RK2 + 2.*RK3 + RK4 ) ;
    	}
    }
    
    for(i=0;i<nbpop;i++) {
      Isyntot[i] = Iext[i] ;
      for(j=0;j<nbpop;j++) 
    	Isyntot[i] += K * J[i][j] * Isyn[i][j] ;
    }
    
    for(i=0;i<nbpop;i++) {      
      r[i] = Func(Isyntot[i]) ;
 
      if(t>=Tst)
    	Mean_Activity[i] += r[i] ; 
    }

    if(t>=Tst) {
      rFile << t ;
      for(i=0;i<nbpop;i++) 
	rFile << " " << r[i] * 1000. ;	
      rFile << endl ;
    }

    if(tw>Tw) {

      cout << " t " << t << " Rates  "  ;
      for(i=0;i<nbpop;i++) {
	cout << Mean_Activity[i] * ( dt / Tw ) * 1000. << " " ;
	Mean_Activity[i] = 0 ;
      }
      cout << endl ;
      cout.flush();
      
      // rFile << t ;
      // for(i=0;i<nbpop;i++) {
      // 	rFile << " " << Mean_Activity[i]*(dt/Tw) ;
      // 	Mean_Activity[i] = 0 ;
      // }
      // rFile << endl ;

      tw=0 ;
    }
    
    // printProgress (percentage) ;
    // updating Time Windows
    if(t>=Tst) 
      tw += dt ;     	  
    
  }//end loop t
  cout << endl ;   
  
  r.clear() ;
  
  Isyn.clear() ;
  Isyntot.clear() ;
  
  rFile.close() ;  
}
