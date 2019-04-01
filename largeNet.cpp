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
  int nbpop, ndI;
  unsigned long N ;
  double K, **J, **Tsyn, g, *Iext, *IextBL, *Crec, Cff = L ;
  
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

  cout << "Membrane time constants : " ;
  for(int i=0;i<nbpop;i++) 
    cout << Tm[i] << " " ;  
  cout << endl ;    
  
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
  // Connectivity Matrix
  ///////////////////////////////////////////////////////////////////

  unsigned long *nbN ;
  nbNeurons(nbpop,N,nbN) ;
  unsigned long *Cpt ; // compteur Cpt0 = Ne, cpt1 = Ne+Ni ...
  cptNeurons(nbpop,nbN,Cpt) ;

  string JpathEE = "../../" ;
  string JpathEI = "../../" ;
  string JpathIE = "../../" ;
  string JpathII = "../../" ;

  string JpathES = "../../" ;
  string JpathIS = "../../" ;
  string JpathSE = "../../" ;
  string JpathSI = "../../" ;

  string JpathSV = "../../" ;
  string JpathVE = "../../" ;
  string JpathVI = "../../" ;
  string JpathVS = "../../" ;

  if(nbpop>=2) {
    Create_Path_Large(JpathEE,nbN[0],K) ;
    Create_Path_Large(JpathEI,nbN[1],K) ;
    Create_Path_Large(JpathIE,nbN[0],K) ;
    Create_Path_Large(JpathII,nbN[1],K) ;
  }
  if(nbpop>=3) {
    Create_Path_Large(JpathES,nbN[0],K) ;
    Create_Path_Large(JpathIS,nbN[1],K) ;
    Create_Path_Large(JpathSE,nbN[0],K) ;
    Create_Path_Large(JpathSI,nbN[1],K) ;
  }
  if(nbpop>=4) {
    Create_Path_Large(JpathSV,nbN[0],K) ;
    Create_Path_Large(JpathVE,nbN[1],K) ;
    Create_Path_Large(JpathVI,nbN[0],K) ;
    Create_Path_Large(JpathVS,nbN[1],K) ;
  }
  if(IF_RING) {
    if(nbpop>=2) {
      CreateDir_SpaceCrec(1,JpathEE,nbN[0],&Crec[0]) ;
      CreateDir_SpaceCrec(1,JpathEI,nbN[1],&Crec[1]) ;
      CreateDir_SpaceCrec(1,JpathIE,nbN[0],&Crec[0]) ;
      CreateDir_SpaceCrec(1,JpathII,nbN[1],&Crec[1]) ;
    }
    if(nbpop>=3) {
      CreateDir_SpaceCrec(1,JpathES,nbN[2],&Crec[2]) ;
      CreateDir_SpaceCrec(1,JpathIS,nbN[2],&Crec[2]) ;
      CreateDir_SpaceCrec(1,JpathSE,nbN[0],&Crec[0]) ;
      CreateDir_SpaceCrec(1,JpathSI,nbN[1],&Crec[1]) ;
    }
    if(nbpop>=4) {
      CreateDir_SpaceCrec(1,JpathSV,nbN[3],&Crec[3]) ;
      CreateDir_SpaceCrec(1,JpathVE,nbN[1],&Crec[0]) ;
      CreateDir_SpaceCrec(1,JpathVI,nbN[0],&Crec[1]) ;
      CreateDir_SpaceCrec(1,JpathVS,nbN[1],&Crec[2]) ;
    }
  }

  int *nbPost_EE, *nbPost_EI , *nbPost_IE, *nbPost_II ;
  unsigned long *idxPost_EE, *idxPost_EI, *idxPost_IE, *idxPost_II ;
  unsigned long *IdPost_EE, *IdPost_EI, *IdPost_IE, *IdPost_II ; 

  int *nbPost_ES, *nbPost_IS , *nbPost_SE, *nbPost_SI ;
  unsigned long *idxPost_ES, *idxPost_IS, *idxPost_SE, *idxPost_SI ;
  unsigned long *IdPost_ES, *IdPost_IS, *IdPost_SE, *IdPost_SI ; 
 
  int *nbPost_SV, *nbPost_VE , *nbPost_VI, *nbPost_VS ;
  unsigned long *idxPost_SV, *idxPost_VE, *idxPost_VI, *idxPost_VS ;
  unsigned long *IdPost_SV, *IdPost_VE, *IdPost_VI, *IdPost_VS ; 

  if(nbpop>=2) {
    Import_Connectivity_Matrix_Large(nbpop,"EE",nbN[0],JpathEE,nbPost_EE,idxPost_EE,IdPost_EE) ;
    Import_Connectivity_Matrix_Large(nbpop,"EI",nbN[1],JpathEI,nbPost_EI,idxPost_EI,IdPost_EI) ;
    Import_Connectivity_Matrix_Large(nbpop,"IE",nbN[0],JpathIE,nbPost_IE,idxPost_IE,IdPost_IE) ;
    Import_Connectivity_Matrix_Large(nbpop,"II",nbN[1],JpathII,nbPost_II,idxPost_II,IdPost_II) ;
  }
  if(nbpop>=3) {
    Import_Connectivity_Matrix_Large(nbpop,"ES",nbN[2],JpathES,nbPost_ES,idxPost_ES,IdPost_ES) ;
    Import_Connectivity_Matrix_Large(nbpop,"IS",nbN[2],JpathIS,nbPost_IS,idxPost_IS,IdPost_IS) ;
    Import_Connectivity_Matrix_Large(nbpop,"SE",nbN[0],JpathSE,nbPost_SE,idxPost_SE,IdPost_SE) ;
    // Import_Connectivity_Matrix_Large(nbpop,"SI",nbN[1],JpathSI,nbPost_SI,idxPost_SI,IdPost_SI) ;
  }
  if(nbpop>=4) {
    Import_Connectivity_Matrix_Large(nbpop,"SV",nbN[0],JpathSV,nbPost_SV,idxPost_SV,IdPost_SV) ;
    Import_Connectivity_Matrix_Large(nbpop,"VE",nbN[1],JpathVE,nbPost_VE,idxPost_VE,IdPost_VE) ;
    Import_Connectivity_Matrix_Large(nbpop,"VI",nbN[0],JpathVI,nbPost_VI,idxPost_VI,IdPost_VI) ;
    Import_Connectivity_Matrix_Large(nbpop,"VS",nbN[1],JpathVS,nbPost_VS,idxPost_VS,IdPost_VS) ;
  }
  
  ///////////////////////////////////////////////////////////////////     
  // Scaling
  //////////////////////////////////////////////////////////////////

  Save_Parameters(dir, nbpop, N, nbN, K, J, Iext, Tm, Tsyn, path) ;

  for(int i=0;i<nbpop;i++) {
    Iext[i] = sqrt(K)*Iext[i]*m0 ;
    IextBL[i] = sqrt(K)*IextBL[i]*m0 ;
    for(int j=0;j<nbpop;j++) 
      J[i][j] = g*J[i][j]/sqrt(K)/Tsyn[i][j] ;
    
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
  if(IF_Prtr)
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
  int i=0,j=0 ;
  unsigned long k=0, l=0 ;
  
  ///////////////////////////////////////////////////////////////////    

  cout << "Initialization" << endl;
  //Using C++11
  random_device rd ;
  default_random_engine gen( rd() ) ;
  uniform_real_distribution<double> unif(0,1) ;

  normal_distribution<double> gaussianV(1,.25);
  normal_distribution<double> gaussianI(0,10.);

  cout << "Check random seed " ;
  for(i=0;i<10;i++)
    cout << unif(gen) << " " ;
  cout << endl ;

  for(i=0;i<nbpop;i++) 
    for(k=0;k<nbN[i];k++) {
      Volt[i][k] = gaussianV(gen) ;
      Isyntot[i][k] = gaussianI(gen) ;
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
	
	// Volt[i][k] = Cst[i] * ( Volt[i][k]-Vr ) + dt * IextBL[i] + Isyntot[i][k] ;
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
	  if(i==0) {
	    for ( l = idxPost_EE[k] ; l < (unsigned long) idxPost_EE[k] + nbPost_EE[k] ; l++) 
	      Isyn[0][i][IdPost_EE[l]] += J[0][i] ; 
	    for ( l = idxPost_IE[k] ; l < (unsigned long) idxPost_IE[k] + nbPost_IE[k] ; l++) 
	      Isyn[1][i][IdPost_IE[l]] += J[1][i] ; 
	    if(nbpop>=3)
	      for ( l = idxPost_SE[k] ; l < (unsigned long) idxPost_SE[k] + nbPost_SE[k] ; l++) 
		Isyn[2][i][IdPost_SE[l]] += J[2][i] ; 
	    if(nbpop>=4)
	      for ( l = idxPost_VE[k] ; l < (unsigned long) idxPost_VE[k] + nbPost_VE[k] ; l++) 
		Isyn[3][i][IdPost_VE[l]] += J[3][i] ; 
	  }

	  if(i==1) {
	    for ( l = idxPost_EI[k] ; l < (unsigned long) idxPost_EI[k] + nbPost_EI[k] ; l++) 
	      Isyn[0][i][IdPost_EI[l]] += J[0][i] ; 
	    for ( l = idxPost_II[k] ; l < (unsigned long) idxPost_II[k] + nbPost_II[k] ; l++) 
	      Isyn[1][i][IdPost_II[l]] += J[1][i] ; 
	    if(nbpop==3)
	      for ( l = idxPost_SI[k] ; l < (unsigned long) idxPost_SI[k] + nbPost_SI[k] ; l++) 
		Isyn[2][i][IdPost_SI[l]] += J[2][i] ; 
	    if(nbpop>=4)
	      for ( l = idxPost_VI[k] ; l < (unsigned long) idxPost_VI[k] + nbPost_VI[k] ; l++) 
		Isyn[3][i][IdPost_VI[l]] += J[3][i] ; 
	  }

	  if(i==2) {
	    for ( l = idxPost_ES[k] ; l < (unsigned long) idxPost_ES[k] + nbPost_ES[k] ; l++) 
	      Isyn[0][i][IdPost_ES[l]] += J[0][i] ; 
	    for ( l = idxPost_IS[k] ; l < (unsigned long) idxPost_IS[k] + nbPost_IS[k] ; l++) 
	      Isyn[1][i][IdPost_IS[l]] += J[1][i] ; 
	    if(nbpop>=4)
	      for ( l = idxPost_VS[k] ; l < (unsigned long) idxPost_VS[k] + nbPost_VS[k] ; l++) 
		Isyn[3][i][IdPost_VS[l]] += J[3][i] ; 
	  }

	  if(i==3) {
	    if(nbpop>=3)
	      for ( l = idxPost_SV[k] ; l < (unsigned long) idxPost_SV[k] + nbPost_SV[k] ; l++) 
		Isyn[2][i][IdPost_SV[l]] += J[2][i] ; 
	  }

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
	    Isyntot[i][k] += Isyn[i][j][k] ; // same as suming over f(t-tspike) see Brette 2006
	    Isyn[i][j][k] = Cstsyn[i][j]*Isyn[i][j][k] ;
	    // Isyntot[i][k] += Isyn[i][j][k] * ( Cst[i] -Cstsyn[i][j] ) * Tsyn[i][j] / (Tm[i] - Tsyn[i][j] ) ;
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

  if(nbpop>=2) {
    delete [] nbPost_EE ;
    delete [] IdPost_EE ;
    delete [] idxPost_EE ;

    delete [] nbPost_EI ;
    delete [] IdPost_EI ;
    delete [] idxPost_EI ;
    
    delete [] nbPost_IE ;
    delete [] IdPost_IE ;
    delete [] idxPost_IE ;
    
    delete [] nbPost_II ;
    delete [] IdPost_II ;
    delete [] idxPost_II ;
  }
  if(nbpop>=3) {

    delete [] nbPost_ES ;
    delete [] IdPost_ES ;
    delete [] idxPost_ES ;
    
    delete [] nbPost_IS ;
    delete [] IdPost_IS ;
    delete [] idxPost_IS ;

    delete [] nbPost_SE ;
    delete [] IdPost_SE ;
    delete [] idxPost_SE ;

  }

  if(nbpop>=4) {

    delete [] nbPost_VS ;
    delete [] IdPost_VS ;
    delete [] idxPost_VS ;
    
    
    delete [] nbPost_VE ;
    delete [] IdPost_VE ;
    delete [] idxPost_VE ;
    
    delete [] nbPost_VI ;
    delete [] IdPost_VI ;
    delete [] idxPost_VI ;
    
    delete [] nbPost_SV ;
    delete [] IdPost_SV ;
    delete [] idxPost_SV ;

  }

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
