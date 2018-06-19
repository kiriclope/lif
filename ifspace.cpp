#include "librairies.h"
#include "Matrix_Utils.h"
#include "Space_Utils.h"
#include "Net_Utils.h"

#define L 2.0*M_PI
#define PULSE 40
#define IF_PULSE false
#define IF_Nk false
#define IF_Iext false
#define PrtrPop 1

#define Vr  0. // Membrane Resting Potential
#define Vth  1. // Voltage Threshold
#define Vpeak 20. // Spikes Peak

#define IF_ADAPTATION 0
#define Tad 100.
#define Gad 1.

clock_t t1=clock();

int main(int argc , char** argv) {

  ///////////////////////////////////////////////////////////////////    
  // parameters 
  ///////////////////////////////////////////////////////////////////    
  
  string dir ;
  int nbpop, N ;
  double K,duration,dt,Tst,Tw,**J,**Tsyn,g,*Iext,*IextBL,*Crec ;
  
  Set_Parameters(argc,argv,dir,nbpop,N,K,duration,dt,Tst,Tw,g,Iext,IextBL) ;

  Crec = new double [nbpop] ;     
  for(int i=0;i<nbpop;i++) 
    Crec[i] = (double) atof(argv[2*nbpop+10+i]) ;

  double Cff = (double) atof(argv[3*nbpop+10]) ; // variance of the wrapped external input
  
  cout << "Spatiality " << endl ;
  cout << "Crec " ;
  for(int i=0;i<nbpop;i++)
    cout << " " << Crec[i] ;
  cout << " Cff " << Cff << endl ;

  ///////////////////////////////////////////////////////////////////    
  int ndI = min(nbpop-1,PrtrPop) ;
  ///////////////////////////////////////////////////////////////////    
  
  cout << "Membrane Time Constantes : " << endl ;
  double *Tm = new double [nbpop]() ;
  for(int i=0;i<nbpop;i++)
    Tm[i] = 10. ;
  
  double *Cst = new double [nbpop]() ;
  if(argv[nbpop+10] != NULL) 
    for(int i=0;i<nbpop;i++) 
      Tm[i] = (double) atof(argv[i+10+nbpop]) ;
  
  for(int i=0;i<nbpop;i++) 
    Cst[i] = exp(-dt/Tm[i]) ;

  for(int i=0;i<nbpop;i++) 
    cout << Tm[i] << " " ;
  cout << endl ;

  double *Cstad = new double [nbpop]() ;
  for(int i=0;i<nbpop;i++) 
    Cstad[i] = exp(-dt/Tad) ;

  ///////////////////////////////////////////////////////////////////    
  
  cout << "Synaptic Time Constantes : " << endl ;
  Import_Synaptic_Parameters(nbpop,Tsyn,dir) ;

  double **Cstsyn = new double *[nbpop]() ;
  for(int i=0;i<nbpop;i++) {
    Cstsyn[i] = new double [nbpop]() ;    
    for(int j=0;j<nbpop;j++)
      Cstsyn[i][j] = exp(-dt/Tsyn[i][j]) ;
  }

  ///////////////////////////////////////////////////////////////////    

  cout << "External Inputs : " ;
  for(int i=0;i<nbpop;i++) 
    cout << Iext[i] << " " ;
  cout << endl ;
  
  cout << "Synaptic Strength : " << endl ;
  Import_Connectivity_Parameters(nbpop,J,dir) ; 
  for(int i=0;i<nbpop;i++) {
    for(int j=0;j<nbpop;j++)
      cout << J[i][j] << " ";
    cout << endl ;
  }

  ///////////////////////////////////////////////////////////////////    
  // Connectivity Matrix
  ///////////////////////////////////////////////////////////////////

  string Jpath = "../../" ;
  Create_Path(nbpop,Jpath,N,K) ;
  CreateDir_SpaceCrec(nbpop,Jpath,N,Crec) ;

  int *nbPostSpace ;
  unsigned long  int *idxPostSpace ;
  int *IdPostSpace ;
  Import_Connectivity_Matrix(nbpop,N,Jpath,nbPostSpace,idxPostSpace,IdPostSpace,IF_Nk) ;

  ///////////////////////////////////////////////////////////////////    
  // Path
  ///////////////////////////////////////////////////////////////////    

  string path = "../" ;
  CreateDir(dir, nbpop, N, K, g, path) ;  
  CreateDir_SpaceCrec(nbpop,path,N,Crec) ;
  CreateDir_SpaceCff(nbpop,path,N,Cff) ;

  if(IF_Iext) {
    cout <<"Perturbed Pop " << ndI << " | Input "<< Iext[ndI]-IextBL[ndI] << endl ;
    CreateDir_Iext(nbpop,ndI,Iext[ndI],path) ;
  }

  ///////////////////////////////////////////////////////////////////    
  // number of neurons
  ///////////////////////////////////////////////////////////////////    

  int* Nk ;
  nbNeurons(nbpop,N,Nk,IF_Nk) ;
    
  int* Cpt ; // compteur Cpt0 = Ne, cpt1 = Ne+Ni ...
  Cpt = new int [nbpop+1]() ;
  
  for(int i=0;i<nbpop+1;i++) { 
    for(int j=0;j<i;j++) 
      Cpt[i] += Nk[j] ; 
    cout <<"Cpt="<< Cpt[i] << " ";
  }
  cout << endl ;

  ///////////////////////////////////////////////////////////////////     
  // Scaling
  //////////////////////////////////////////////////////////////////

  duration = duration*1000. ;
  Tst = Tst*1000. ;
  Tw = Tw*1000. ;

  Save_Parameters(dir, nbpop, N, Nk, K, J, duration, dt, Tst, Tw, Iext, Tm, Tsyn, path) ;

  for(int i=0;i<nbpop;i++) {
    Iext[i] = sqrt(K)*Iext[i]*m0 ;
    IextBL[i] = sqrt(K)*IextBL[i]*m0 ;
    for(int j=0;j<nbpop;j++) {
      J[i][j] = g*J[i][j]/sqrt(K)/Tsyn[i][j] ;
    }
  }

  ///////////////////////////////////////////////////////////////////    
  // Dynamics of the network : I&F
  ///////////////////////////////////////////////////////////////////    
  
  vector<double> Mean_Activity(nbpop) ;
  vector<vector<double> > Idv_Inputs(nbpop,vector<double>(N)) ;
  vector<vector<double> > Idv_Activity(nbpop,vector<double>(N)) ;

  vector<vector<double> > Jext(nbpop,vector<double>(N)) ;

  vector<vector<vector<double> > >Isyn(nbpop,vector<vector<double> >(nbpop)) ; 

  vector<vector<double> > Isyntot(nbpop,vector<double>(N)) ; 
  vector<vector<double> > Volt(nbpop,vector<double>(N)) ; 

  vector<vector<double> > Iad(nbpop,vector<double>(N)) ;

  for (int i=0;i<nbpop;i++) {
    Idv_Activity[i].resize(Nk[i]) ;
    Volt[i].resize(Nk[i]) ;
    Isyntot[i].resize(Nk[i]) ;
    Iad[i].resize(Nk[i]) ;
    Jext[i].resize(Nk[i]) ;
    for (int j=0;j<nbpop;j++) 
      Isyn[i][j].resize(Nk[i]) ;    
  }

  External_Input(nbpop,N,Nk,K,Cff,Iext,IextBL,ndI,Jext,L,path) ;


  ///////////////////////////////////////////////////////////////////    

  string strMean = path + "/Mean.txt" ;
  ofstream fMean(strMean.c_str(), ios::out | ios::ate);

  string strInput = path + "/Input.txt" ;
  ofstream fInput(strInput.c_str(), ios::out | ios::ate);

  string strVoltage = path + "/Voltage.txt" ;
  ofstream fVoltage(strVoltage.c_str(), ios::out | ios::ate);

  string strRaster = path + "/Raster.txt" ;
  ofstream fRaster(strRaster.c_str(), ios::out | ios::ate);

  string strIdvInputs = path + "/IdvInputs.txt" ;
  ofstream fIdvInputs(strIdvInputs.c_str(), ios::out | ios::ate);

  string strIdvRates = path + "/IdvRates.txt" ;
  ofstream fIdvRates(strIdvRates.c_str(), ios::out | ios::ate);

  ///////////////////////////////////////////////////////////////////    

  double tw = 0. ; //Time window
  int i=0,j=0,k=0,m=0;
  unsigned long int l=0;
  
  ///////////////////////////////////////////////////////////////////    

  cout << "Initialization" << endl;
  //Using C++11
  random_device rd ;
  default_random_engine gen( rd() ) ;
  uniform_real_distribution<double> unif(0,1) ;
  normal_distribution<double> gaussian(1.,.25);

  cout << "Check random seed " ;
  for(int i=0;i<10;i++)
    cout << unif(gen) << " " ;
  cout << endl ;

  for(i=0;i<nbpop;i++) 
    for(k=0;k<Nk[i];k++) {
      Volt[i][k] = gaussian(gen) ;
      Isyntot[i][k] = gaussian(gen) ;
    }
  
  ///////////////////////////////////////////////////////////////////    
  // Main Loop
  ///////////////////////////////////////////////////////////////////     

  double RK1=0,RK2=0,RK3=0,RK4=0 ;
  double Vold=0, percentage=0 ;
  double Tl = 2000. ;

  cout << "Main loop :" ;
  cout << " duration " << duration << " | dt " << dt ;
  cout << " | Tst " << Tst << " | Tw " << Tw << " | Tl " << Tl << endl ;

  for (double t=0.; t<=duration; t+=dt) {  

    percentage = t/duration ;
    
    if(t>=Tst && t<=Tst+Tl ) fVoltage << t-Tst ; // Adding Spikes by hand 
    
    for(i=0;i<nbpop;i++) {    
      for (k=0; k<Nk[i]; k++) {
	
	// if(t>=Tst && k<10) Cvl_Idv_Activity[i][k] = Cvl_Idv_Activity[i][k]*expalpha ;

	//Updating Voltage 
	
	///////////////////////////////////////////////////////////////////    
	// Standard Euler Algorithm
	///////////////////////////////////////////////////////////////////    
	
	// Volt[i][k] = Cst[i]*(Volt[i][k]-Vr)+dt*Isyntot[i][k] ;
	// Volt[i][k]=(1.-dt/Tm[i])*Volt[i][k]+dt*(Isyntot[i][k]+Vr/Tm[i]);
	
	///////////////////////////////////////////////////////////////////    
	// RK2
	///////////////////////////////////////////////////////////////////    

	// RK1 = -(Volt[i][k]-Vr)/Tm[i] + Isyntot[i][k] ;
	// RK2 = -(Volt[i][k]-Vr+dt*RK1)/Tm[i] + Isyntot[i][k] ;
	// Volt[i][k] = Volt[i][k] + dt/2.*(RK1 + RK2 );
	
	///////////////////////////////////////////////////////////////////    
	// RK4
	///////////////////////////////////////////////////////////////////    

	// Vold = Volt[i][k] ;

	RK1 = -(Volt[i][k]-Vr)/Tm[i] + Isyntot[i][k] ;
	RK2 = -(Volt[i][k]-Vr+dt/2.*RK1)/Tm[i] + Isyntot[i][k] ;
	RK3 = -(Volt[i][k]-Vr+dt/2.*RK2)/Tm[i] + Isyntot[i][k] ;
	RK4 = -(Volt[i][k]-Vr+dt*RK3)/Tm[i] + Isyntot[i][k] ;
	Volt[i][k] = Volt[i][k] + dt/6.*(RK1 + 2.*RK2 + 2.*RK3 + RK4) ;
		
	///////////////////////////////////////////////////////////////////    
	// IF Spike
	///////////////////////////////////////////////////////////////////    

	//if Spike reset V and update postsynaptic currents
	if (Volt[i][k]>=Vth) {
	  Volt[i][k] = Vr ;
	  
	  if(IF_ADAPTATION & i==0)
	    Iad[i][k] += Gad/Tad ;

	  // Improved accuracy
	  // Volt[i][k] = (Volt[i][k]-Vth)*(1.+dt/Tm[i]*(Vold-Vr)/(Volt[i][k]-Vold)) + Vr ;
	  
	  if(t>=Tst) {
	    Mean_Activity[i] += 1. ;
	    Idv_Inputs[i][k] += Isyntot[i][k] ;
	    Idv_Activity[i][k] += 1. ;
	    if(k<10 && t<=Tst+Tl) fVoltage << " " << Vpeak ;
	  }

	  // Updating Postsynaptics inputs 
	  // for(j=0;j<nbpop;j++) // i Pres to j Post => Jji 
	  //   if(J[j][i]!=0)
	  //     for (l=idxPost[k+Cpt[i]]; l<idxPost[k+Cpt[i]]+nbPost[k+Cpt[i]]; l++) 
	  // 	if(IdPost[l]>=Cpt[j] && IdPost[l]<Cpt[j+1]) 
	  // 	  Isyn[j][i][IdPost[l]-Cpt[j]] += J[j][i] ; 

	  for(j=0;j<nbpop;j++) // i Pres to j Post => Jji 
	    if(J[j][i]!=0)
	      for (l=idxPostSpace[k+Cpt[i]]; l<idxPostSpace[k+Cpt[i]]+nbPostSpace[k+Cpt[i]]; l++) 
		if(IdPostSpace[l]>=Cpt[j] && IdPostSpace[l]<Cpt[j+1]) 
		  Isyn[j][i][IdPostSpace[l]-Cpt[j]] += J[j][i] ; 
	  
	  //Updating ISI  
	  if(t>=Tst && t<=Tst+Tl)
	    fRaster << k+Cpt[i] << " " << t-Tst << endl ;
	  
	} // endif spike 
	else
	  if(t>=Tst and k<10 && t<=Tst+Tl) fVoltage << " " << Volt[i][k] ;
      } // endfor neurons
    } // endfor populations
    
    if(t>=Tst && t<=Tst+Tl) fVoltage << endl ;
    
    // Updating Total Synaptic Input to each neurons in each population 
    for(i=0;i<nbpop;i++) 
      for (k=0;k<Nk[i];k++) 
	Isyntot[i][k] = IextBL[i] + Jext[i][k] ;

    for(i=0;i<nbpop;i++) 
      for(j=0;j<nbpop;j++) 
	if(J[i][j]!=0) 
	  for (k=0;k<Nk[i];k++) {
	    Isyn[i][j][k] = Cstsyn[i][j]*Isyn[i][j][k] ; // same as suming over f(t-tspike) see Brette 2006
	    Isyntot[i][k] += Isyn[i][j][k] ;
	    // Isyn[i][j][k] = 0 ; // instantaneous synapse 
	  }
    
    // Adaptation
    if(IF_ADAPTATION)
      for(i=0;i<nbpop;i++) 
	if(i==0)
	  for (k=0;k<Nk[i];k++) {
	    Iad[i][k] = Cstad[i]*Iad[i][k] ; 
	    Isyntot[i][k] -= Iad[i][k] ;
	  }

    // Writing to file  
    if(t>=Tst && t<=Tst+Tl) {
      
      fInput << t-Tst ;
      for(int i=0;i<nbpop;i++) {
      	fInput << " " << (Iext[i]+Isyn[i][0][0])*Tm[i] ;
      	for(int j=1;j<nbpop;j++) 
      	  fInput << " " << Isyn[i][j][0]*Tm[i] ;
      }
      fInput << endl;  
      
    } //endif 
    
    if(tw>=Tw) {

      cout << " t " << t-Tst << " Rates " ;
      for(i=0;i<nbpop;i++) {
    	cout << " " << Mean_Activity[i]/tw/(double)Nk[i]*1000. ;
	Mean_Activity[i] = 0 ;
      }
      // cout << endl ;

      cout.flush();
      
      if(t<=Tst+Tl) {

	fMean << t-Tst ;
	for(i=0;i<nbpop;i++) {
	  fMean << " " << Mean_Activity[i]/tw/(double)Nk[i]*1000. ;
	  Mean_Activity[i] = 0 ;
	}
	fMean << endl ;
	
	fIdvInputs << t-Tst ;
	for (int i=0 ;i<nbpop;i++) 
	  for (int k=0 ;k<Nk[i];k++) {
	    fIdvInputs <<" "<< Idv_Inputs[i][k]/tw*1000. ;
	    Idv_Inputs[i][k] = 0 ;
	  }
	fIdvInputs << endl ;
      }

      fIdvRates << t-Tst ;
      for (int i=0 ;i<nbpop;i++) 
	for (int k=0 ;k<Nk[i];k++) {
	  fIdvRates <<" "<< Idv_Activity[i][k]/tw*1000. ;
	  Idv_Activity[i][k] = 0 ;
	}
      fIdvRates << endl ;

      // if(IF_OPSIN)
      // 	External_Input(nbpop,N,Nk,K,Cff,Iext,IextBL,ndI,Jext,L,path) ;

      tw=0;
      
    }//ENDIF
    
    printProgress (percentage) ;

    // updating Time Windows
    if(t>=Tst)
      tw += dt ;

  } //ENDMAINLOOP
  
  cout << endl ;

  ///////////////////////////////////////////////////////////////////    

  // delete[] nbPost ;
  // delete[] IdPost ;
  // delete[] idxPost ;

  delete[] nbPostSpace ;
  delete[] IdPostSpace ;
  delete[] idxPostSpace ;

  delete [] Nk ;
  delete [] Cpt ;

  delete [] Iext ;
  delete [] J ;

  delete [] Tm ;
  delete [] Tsyn ;

  Isyn.clear() ;
  Isyntot.clear() ;
  Volt.clear() ;

  Jext.clear() ;
  Mean_Activity.clear() ;
  Idv_Inputs.clear() ;
  Idv_Activity.clear() ;

  ///////////////////////////////////////////////////////////////////    

  fMean.close();
  fInput.close();
  fVoltage.close();
  fRaster.close();
  fIdvInputs.close();
  fIdvRates.close();

  ///////////////////////////////////////////////////////////////////    

  cout << "Simulation Done !" << endl ;

  clock_t t2=clock();  
  int HOURS=0,MIN=0,SEC=0;
  string str_TIME = path + "/CPU_TIME.txt" ; 
  ofstream TIME(str_TIME.c_str(), ios::out | ios::ate);

  SEC = (t2-t1)/CLOCKS_PER_SEC ;
  HOURS = SEC/3600 ;
  MIN = SEC/60 ;
  SEC = SEC % 60 ;
  cout << "Elapsed Time = " << HOURS << "h " << MIN << "m " << SEC << "s" << endl;
  TIME << "Elapsed Time = " << HOURS << "h " << MIN << "m " << SEC << "s" << endl;
  TIME.close() ;
  return 0;

}
