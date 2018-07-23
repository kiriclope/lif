#include "librairies.h"
#include "GlobalVars.h"
#include "Net_Utils.h"
#include "Matrix_Utils.h"
#include "Space_Utils.h"

clock_t t1=clock();

int main(int argc , char** argv) {

  ///////////////////////////////////////////////////////////////////    
  // parameters 
  ///////////////////////////////////////////////////////////////////    
  
  string dir ;
  int nbpop, N, ndI;
  double K, **J, **Tsyn, g, *Iext, *IextBL, *Crec, Cff ;
  
  Set_Parameters(argc,argv,dir,nbpop,N,K,g,Iext,IextBL) ; 

  if(IF_Prtr) {
    ndI = min(nbpop-1,PrtrPop) ;
    Iext[ndI] += (double) atof(argv[6]) ;
  }

  if(IF_RING)
    getCrecCff(argv,nbpop,Crec,Cff) ;
  
  cout << "Synaptic Strength : " << endl ;
  Import_Connectivity_Parameters(nbpop,J,dir) ; 
  for(int i=0;i<nbpop;i++) {
    for(int j=0;j<nbpop;j++)
      cout << J[i][j] << " ";
    cout << endl ;
  }

  ///////////////////////////////////////////////////////////////////    
  // Time Csts
  ///////////////////////////////////////////////////////////////////    
  
  double *Cst = new double [nbpop]() ;
  for(int i=0;i<nbpop;i++) 
    Cst[i] = exp(-dt/Tm[i]) ;
   
  Import_Synaptic_Parameters(nbpop,Tsyn,dir) ;

  double **Cstsyn = new double *[nbpop]() ;
  for(int i=0;i<nbpop;i++) {
    Cstsyn[i] = new double [nbpop]() ;    
    for(int j=0;j<nbpop;j++)
      Cstsyn[i][j] = exp(-dt/Tsyn[i][j]) ;
  }

  ///////////////////////////////////////////////////////////////////    
  // Connectivity Matrix
  ///////////////////////////////////////////////////////////////////

  string Jpath = "../../" ;
  Create_Path(nbpop,Jpath,N,K) ;
  if(IF_RING)
    CreateDir_SpaceCrec(nbpop,Jpath,N,Crec) ;

  int *nbPost ;
  unsigned long  int *idxPost ;
  int *IdPost ;
  Import_Connectivity_Matrix(nbpop,N,Jpath,nbPost,idxPost,IdPost) ;

  ///////////////////////////////////////////////////////////////////    
  // Path
  ///////////////////////////////////////////////////////////////////    

  string path = "../" ;
  CreateDir(dir, nbpop, N, K, g, path) ; 

  if(IF_RING) {
    CreateDir_SpaceCrec(nbpop,path,N,Crec) ;
    CreateDir_SpaceCff(nbpop,path,N,Cff) ;
  }

  string popList[4] = {"E","I","S","V"} ;  
  if(nbpop==1) 
    popList[0] = "I" ;

  if(IF_Prtr) {
    cout <<"Perturbed Pop " << popList[ndI] << " | Input "<< Iext[ndI]-IextBL[ndI] << endl ;
    CreateDir_Iext(nbpop,ndI,Iext[ndI],path) ;
  }

  ///////////////////////////////////////////////////////////////////    
  // number of neurons
  ///////////////////////////////////////////////////////////////////    

  int* nbN ;
  nbNeurons(nbpop,N,nbN) ;
  int* Cpt ; // compteur Cpt0 = Ne, cpt1 = Ne+Ni ...
  cptNeurons(nbpop,nbN,Cpt) ;

  ///////////////////////////////////////////////////////////////////     
  // Scaling
  //////////////////////////////////////////////////////////////////

  Save_Parameters(dir, nbpop, N, nbN, K, J, Iext, Tm, Tsyn, path) ;

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

  for (int i=0;i<nbpop;i++) {
    Idv_Activity[i].resize(nbN[i]) ;
    Volt[i].resize(nbN[i]) ;
    Isyntot[i].resize(nbN[i]) ;
    Jext[i].resize(nbN[i]) ;
    for (int j=0;j<nbpop;j++) 
      Isyn[i][j].resize(nbN[i]) ;    
  }
  
  // if(IF_Prtr & !IF_RING)
  //   External_Input(nbpop,N,nbN,K,100,Iext,IextBL,ndI,Jext,path) ;
  if(IF_Prtr & IF_RING)
    External_Input(nbpop,N,nbN,K,Cff,Iext,IextBL,ndI,Jext,path) ;

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
    for(k=0;k<nbN[i];k++) {
      Volt[i][k] = gaussian(gen) ;
      Isyntot[i][k] = gaussian(gen) ;
    }
  
  ///////////////////////////////////////////////////////////////////    
  // Main Loop
  ///////////////////////////////////////////////////////////////////     

  double RK1=0,RK2=0,RK3=0,RK4=0 ;
  double Vold=0, percentage=0 ;

  cout << "Main loop :" ;
  cout << " duration " << duration << " | dt " << dt ;
  cout << " | Tst " << Tst << " | Tw " << Tw << " | Tl " << Tl << endl ;

  for (double t=0.; t<=duration; t+=dt) {  

    percentage = t/duration ;
    
    if(t>=Tst && t<=Tst+Tl ) fVoltage << t-Tst ; // Adding Spikes by hand 
    
    for(i=0;i<nbpop;i++) {    
      for (k=0; k<nbN[i]; k++) {
	
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
	  
	  // Improved accuracy
	  // Volt[i][k] = (Volt[i][k]-Vth)*(1.+dt/Tm[i]*(Vold-Vr)/(Volt[i][k]-Vold)) + Vr ;
	  
	  if(t>=Tst) {
	    Mean_Activity[i] += 1. ;
	    Idv_Inputs[i][k] += Isyntot[i][k] ;
	    Idv_Activity[i][k] += 1. ;
	    if(k<10 && t<=Tst+Tl) fVoltage << " " << Vpeak ;
	  }

	  // Updating Postsynaptics inputs 
	  for(j=0;j<nbpop;j++) // i Pres to j Post => Jji 
	    if(J[j][i]!=0)
	      for (l=idxPost[k+Cpt[i]]; l<idxPost[k+Cpt[i]]+nbPost[k+Cpt[i]]; l++) 
		if(IdPost[l]>=Cpt[j] && IdPost[l]<Cpt[j+1]) 
		  Isyn[j][i][IdPost[l]-Cpt[j]] += J[j][i] ; 
	  
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
      for (k=0;k<nbN[i];k++) 
	Isyntot[i][k] = IextBL[i] + Jext[i][k] ;

    for(i=0;i<nbpop;i++) 
      for(j=0;j<nbpop;j++) 
	if(J[i][j]!=0) 
	  for (k=0;k<nbN[i];k++) {
	    Isyn[i][j][k] = Cstsyn[i][j]*Isyn[i][j][k] ; // same as suming over f(t-tspike) see Brette 2006
	    Isyntot[i][k] += Isyn[i][j][k] ;
	    // Isyn[i][j][k] = 0 ; // instantaneous synapse 
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
    	cout << " " << Mean_Activity[i]/tw/(double)nbN[i]*1000. ;
	Mean_Activity[i] = 0 ;
      }
      // cout << endl ;

      cout.flush();
      
      if(t<=Tst+Tl) {

	fMean << t-Tst ;
	for(i=0;i<nbpop;i++) {
	  fMean << " " << Mean_Activity[i]/tw/(double)nbN[i]*1000. ;
	  Mean_Activity[i] = 0 ;
	}
	fMean << endl ;
	
	fIdvInputs << t-Tst ;
	for (int i=0 ;i<nbpop;i++) 
	  for (int k=0 ;k<nbN[i];k++) {
	    fIdvInputs <<" "<< Idv_Inputs[i][k]/tw*1000. ;
	    Idv_Inputs[i][k] = 0 ;
	  }
	fIdvInputs << endl ;
      }

      fIdvRates << t-Tst ;
      for (int i=0 ;i<nbpop;i++) 
	for (int k=0 ;k<nbN[i];k++) {
	  fIdvRates <<" "<< Idv_Activity[i][k]/tw*1000. ;
	  Idv_Activity[i][k] = 0 ;
	}
      fIdvRates << endl ;

      // if(IF_OPSIN)
      // 	External_Input(nbpop,N,nbN,K,Cff,Iext,IextBL,ndI,Jext,L,path) ;

      tw=0;
      
    }//ENDIF
    
    printProgress (percentage) ;

    // updating Time Windows
    if(t>=Tst)
      tw += dt ;

  } //ENDMAINLOOP
  
  cout << endl ;

  ///////////////////////////////////////////////////////////////////    

  delete[] nbPost ;
  delete[] IdPost ;
  delete[] idxPost ;

  delete [] nbN ;
  delete [] Cpt ;

  delete [] Iext ;
  delete [] IextBL ;
  delete [] J ;
 
  delete [] Crec ;

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
