#include "librairies.h"
#include "GlobalVars.h" 
#include "Net_Utils.h"
#include "Matrix_Utils.h"
#include "Space_Utils.h"
#include "MFrates.h"
clock_t t1=clock();

int main(int argc , char** argv) {

  ///////////////////////////////////////////////////////////////////    
  // parameters 
  ///////////////////////////////////////////////////////////////////    
  
  string dir ;
  int nbpop, ndI ; 
  unsigned long N ;  
  double dt, K, **J, **Tsyn, g, *Iext, *IextBL, *Crec, Cff = 4.0*L ; 
  double *RatesMF ;

  dt = DT ;
  if(IF_BENCHMARK)
    dt = (double) atof(argv[6]) ;
  
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

  MF_Rates(nbpop,IextBL,J,RatesMF) ;

  cout << "MF Rates: " ;
  for(int i=0;i<nbpop;i++) 
    cout << RatesMF[i] << " ";
  cout << endl ;

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
    for(int j=0;j<nbpop;j++) { 
      Tsyn[i][j] = Tsyn[i][j] ; 
      Cstsyn[i][j] = exp(-dt/Tsyn[i][j]) ; 
    }
  }

  ///////////////////////////////////////////////////////////////////    
  // Connectivity Matrix
  ///////////////////////////////////////////////////////////////////

  string Jpath = "../../" ;
  Create_Path(nbpop,Jpath,N,K) ;
  if(IF_RING) 
    CreateDir_SpaceCrec(nbpop,Jpath,N,Crec) ;

  int *nbPost ;
  unsigned long *idxPost ;
  unsigned long *IdPost ; 
  Import_Connectivity_Matrix(nbpop,N,Jpath,nbPost,idxPost,IdPost) ; 
  
  ///////////////////////////////////////////////////////////////////
  // Path 
  ///////////////////////////////////////////////////////////////////    

  string path = "../" ; 
  CreateDir(dir, nbpop, N, K, g, path) ; 

  if(IF_SHARED)
    CreateDirShared(path) ;

  if(IF_RING) { 
    CreateDir_SpaceCrec(nbpop,path,N,Crec) ; 
    if(IF_Prtr)
      CreateDir_SpaceCff(nbpop,path,N,Cff) ; 
  } 

  string popList[4] = {"E","I","S","V"} ; 
  if(nbpop==1) 
    popList[0] = "I" ; 

  if(IF_Prtr) {
    cout <<"Perturbed Pop " << popList[ndI] << " | Input "<< Iext[ndI]-IextBL[ndI] << endl ;
    CreateDir_Iext(nbpop,ndI,Iext[ndI],path) ;
  }

  if(IF_JabLoop) { 
    J[Ax][By] = (double) atof(argv[ IF_RING * (nbpop + 1) + 6 + IF_Prtr ]) ; 
    cout <<" J_ " << popList[Ax] << popList[By] << " " << J[Ax][By] << endl ;
    CreateDir_JabLoop(path, J[Ax][By]) ;
  }

  ///////////////////////////////////////////////////////////////////    
  // number of neurons
  ///////////////////////////////////////////////////////////////////    

  unsigned long *nbN ;
  nbNeurons(nbpop,N,nbN) ;
  unsigned long *Cpt ; // compteur Cpt0 = Ne, cpt1 = Ne+Ni ...
  cptNeurons(nbpop,nbN,Cpt) ;

  ///////////////////////////////////////////////////////////////////     
  // Scaling
  //////////////////////////////////////////////////////////////////

  Save_Parameters(dir, nbpop, N, nbN, dt, K, J, Iext, Tm, Tsyn, path) ;

  for(int i=0;i<nbpop;i++) {
    if(~IF_TRANSIENT_IEXT) 
      Iext[i] = sqrt(K) * Iext[i] * m0 * (Vth-Vr) ; 
    IextBL[i] = sqrt(K) * IextBL[i] * m0 * (Vth-Vr) ; 
    for(int j=0;j<nbpop;j++) 
      J[i][j] = g * J[i][j] / sqrt(K) / Tsyn[i][j] * (Vth-Vr) ; 
  } 
  

  ///////////////////////////////////////////////////////////////////    
  // Dynamics of the network : I&F
  ///////////////////////////////////////////////////////////////////    
  
  vector<double> Mean_Activity(nbpop) ; 
  vector<vector<double> > Idv_Inputs(nbpop,vector<double>(N)) ; 
  vector<vector<double> > Idv_Activity(nbpop,vector<double>(N)) ; 

  vector<vector<double> > Jext(nbpop,vector<double>(N)) ; 

  vector<vector<vector<double> > >Isyn(nbpop,vector<vector<double> >(nbpop)) ; 
  vector<vector<vector<double> > >Mean_Isyn(nbpop,vector<vector<double> >(nbpop)) ; 

  vector<vector<double> > Isyntot(nbpop,vector<double>(N)) ; 
  vector<vector<double> > IsyntotRK2(nbpop,vector<double>(N)) ; 

  vector<vector<double> > Volt(nbpop,vector<double>(N)) ; 

  vector<vector<double> > Vold(nbpop,vector<double>(N)) ; 
  vector<vector<double> > tspk(nbpop,vector<double>(N)) ; 

  for (int i=0;i<nbpop;i++) {
    Idv_Activity[i].resize(nbN[i]) ;
    Volt[i].resize(nbN[i]) ;

    Vold[i].resize(nbN[i]) ;
    tspk[i].resize(nbN[i]) ;

    Isyntot[i].resize(nbN[i]) ;
    IsyntotRK2[i].resize(nbN[i]) ;

    Jext[i].resize(nbN[i]) ;
    for (int j=0;j<nbpop;j++) {
      Isyn[i][j].resize(nbN[i]) ;
      Mean_Isyn[i][j].resize(nbN[i]) ;       
    }
  }
  
  if(IF_Prtr)
    if(DIM==1) 
      External_Input(nbpop,N,nbN,K,Cff,Iext,IextBL,ndI,Jext,path) ;
    else
      External_Input2D(nbpop,N,nbN,K,Cff,Iext,IextBL,ndI,Jext,path) ;

  ///////////////////////////////////////////////////////////////////    

  if(IF_BENCHMARK) {

    char cdt[10] ;
    sprintf(cdt,"%0.3f",dt) ; 
    string sdt = string(cdt) ; 
    
    path += "/BenchMark/dt" + sdt ; 

    string mkdirp = "mkdir -p " ; 
    mkdirp += path ; 

    const char * cmd = mkdirp.c_str() ; 
    const int dir_err = system(cmd) ; 
    
    if(-1 == dir_err) 
      cout << "error creating directories" << endl ; 
    cout << path << endl ; 
    
  }

  string strMean = path + "/Mean.txt" ;
  ofstream fMean(strMean.c_str(), ios::out | ios::ate);
  
  string strMeanInputs = path + "/MeanInputs.txt" ;
  ofstream fMeanInputs(strMeanInputs.c_str(), ios::out | ios::ate);

  // string strInput = path + "/Input.txt" ;
  // ofstream fInput(strInput.c_str(), ios::out | ios::ate);

  string strVoltage = path + "/Voltage.txt" ;
  ofstream fVoltage(strVoltage.c_str(), ios::out | ios::ate);

  string strRaster = path + "/Raster.txt" ;
  ofstream fRaster(strRaster.c_str(), ios::out | ios::ate);

  string strIdvInputs = path + "/IdvInputs.txt" ;
  ofstream fIdvInputs(strIdvInputs.c_str(), ios::out | ios::ate);

  string strIdvRates = path + "/IdvRates.txt" ;
  ofstream fIdvRates(strIdvRates.c_str(), ios::out | ios::ate);

  ///////////////////////////////////////////////////////////////////    

  double tw=0., tc=0.  ; //Time window
  int i=0,j=0 ;
  unsigned long k=0, l=0 ;
  int IF_TRANSIENT = 0 ;
  
  ///////////////////////////////////////////////////////////////////    

  cout << "Initialization" << endl;
  //Using C++11
  random_device rd ;
  default_random_engine gen( rd() ) ;
  uniform_real_distribution<double> unif(0,1) ;

  uniform_real_distribution<double> unifV(Vr,Vth) ;
  uniform_real_distribution<double> unifI(Vr,Vth) ;

  normal_distribution<double> gaussianV( (Vth-Vr)/2.0, (Vth-Vr)/8.0 ) ; 
  normal_distribution<double> gaussianI( 0, (Vth-Vr)/8.0 ) ; 

  cout << "Check random seed " ; 
  for(i=0;i<10;i++) 
    cout << unif(gen) << " " ; 
  cout << endl ; 
  
  double I0 = 1.0 ;

  for(i=0;i<nbpop;i++) 
    for(k=0;k<nbN[i];k++) { 
      // Volt[i][k] = unifV(gen) ; 
      // Isyntot[i][k] = unifI(gen) ; 
      // Volt[i][k] = RatesMF[i] + gaussianV(gen) ; 
      // Volt[i][k] = gaussianV(gen) ; 

      // Volt[i][k] = Vr + Tm[i] * I0 * ( 1.0 - exp( ( (double) i-1.0 ) / (double) nbN[i] * log(1.0 - (Vth - Vr) / Tm[i] / I0 ) ) ) ; 
      Volt[i][k] = Tm[i]*Vr ; 
      
      for(j=0;j<nbpop;j++) { 
	// Isyn[i][j][k] = gaussianV(gen) ; 
	Isyn[i][j][k] = ( K * J[i][j] * Tsyn[i][j] * RatesMF[j] + gaussianI(gen) / Tm[i] ) / (double) nbN[i] ; 
      }
    }
  
  cout << "Scheme: " ;
  if(IF_EULER) cout << "EULER " ;
  if(IF_RK2) cout << "RK2 " ;
  if(IF_INTERPOLATION) cout << "with interpolation" ; 
  cout << endl ;
  
  /////////////////////////////////////////////////////////////////// 
  // Main Loop 
  ///////////////////////////////////////////////////////////////////     

  double RK1=0,RK2=0,RK3=0,RK4=0 ;
  double percentage=0 ;
  double tIext=0 ;
  
  cout << "Main loop :" ;
  cout << " duration " << duration << " | dt " << dt ;
  cout << " | Tst " << Tst << " | Tw " << Tw << " | Tl " << Tl << endl ;

  for (double t=0.; t<=duration+Tst+dt; t+=dt) { 

    percentage = t/(duration+Tst) ;
    
    if(t>=duration+Tst-Tl ) fVoltage << t-duration-Tst+Tl ; // Adding Spikes by hand 
    
    // Updating Total Synaptic Input to each neurons in each population 
    for(i=0;i<nbpop;i++) 
      for (k=0;k<nbN[i];k++) {
	Isyntot[i][k] = Iext[i] ; 
	IsyntotRK2[i][k] = Iext[i] ; 
	// Isyntot[i][k] = Iext[i] + Jext[i][k] ; 
    	// Isyntot[i][k] = max(Iext[i] + Jext[i][k] , 0.0 ) ;
	// Isyntot[i][k] = IextBL[i] + Jext[i][k] + .5*( 1.0 + cos(2*M_PI*.003*t) ) ;
      }
    
    for(i=0;i<nbpop;i++) 
      for(j=0;j<nbpop;j++) 
    	if(J[i][j]!=0) 
    	  for (k=0;k<nbN[i];k++) {
    	    // Mean_Isyn[i][j][k] += Isyn[i][j][k] ; 
	    Isyn[i][j][k] = Cstsyn[i][j]*Isyn[i][j][k] ; 
    	    Isyntot[i][k] += Isyn[i][j][k] ; 
	    IsyntotRK2[i][k] += Cstsyn[i][j]*Isyn[i][j][k] ; // Total current at t+dt as if no spike were emitted 
    	    // Isyn[i][j][k] = 0 ; // instantaneous synapse 
    	  }
    
    //Updating Voltage    
    for(i=0;i<nbpop;i++) { 
      for (k=0; k<nbN[i]; k++) { 
	Vold[i][k] = Volt[i][k] ; 
		
	// if(t>=Tst && k<10) Cvl_Idv_Activity[i][k] = Cvl_Idv_Activity[i][k]*expalpha ;	
	  
	///////////////////////////////////////////////////////////////////    
	// Standard Euler Algorithm 
	///////////////////////////////////////////////////////////////////    	
	if(IF_EULER){
	  Volt[i][k] = ( 1.0 - dt / Tm[i] ) * Volt[i][k] + dt * ( Isyntot[i][k] + Vr / Tm[i] ) ; 
	  // Volt[i][k] = Cst[i] * ( Volt[i][k]-Vr ) + dt * Isyntot[i][k] ;	  
	}
	else{
	  ///////////////////////////////////////////////////////////////////    
	  // RK2
	  ///////////////////////////////////////////////////////////////////    
	  if(IF_RK2) {
	    RK1 = -(Volt[i][k]-Vr)/Tm[i] + Isyntot[i][k] ; 
	    RK2 = -(Volt[i][k]-Vr+dt*RK1)/Tm[i] + IsyntotRK2[i][k] ; 
	    Volt[i][k] = Volt[i][k] + dt/2.0 * ( RK1 + RK2 ) ; 	
	  }
	}
	  
	///////////////////////////////////////////////////////////////////    
	// RK4
	///////////////////////////////////////////////////////////////////    

	// RK1 = -(Volt[i][k]-Vr)/Tm[i] + Isyntot[i][k] ; 
	// RK2 = -(Volt[i][k]-Vr+dt/2.*RK1)/Tm[i] + Isyntot[i][k] ; 
	// RK3 = -(Volt[i][k]-Vr+dt/2.*RK2)/Tm[i] + Isyntot[i][k] ; 
	// RK4 = -(Volt[i][k]-Vr+dt*RK3)/Tm[i] + Isyntot[i][k] ; 
	// Volt[i][k] = Volt[i][k] + dt/6.*(RK1 + 2.*RK2 + 2.*RK3 + RK4) ; 
		
	///////////////////////////////////////////////////////////////////    
	// IF Spike
	///////////////////////////////////////////////////////////////////    

	//if Spike reset V and update postsynaptic currents
	if (Volt[i][k]>=Vth) {
	  
	  ///////////////////////////////////////////////////////////////////    
	  // Improved accuracy
	  /////////////////////////////////////////////////////////////////// 
	  if(IF_INTERPOLATION) {
	    tspk[i][k] = t + dt * (Vth-Vold[i][k]) / (Volt[i][k] - Vold[i][k]) ; 
	    Volt[i][k] = ( Volt[i][k] - Vth ) * ( 1.0 + dt/Tm[i] * ( Vold[i][k] - Vr) / ( Volt[i][k] - Vold[i][k] ) ) + Vr ; 
	  } 
	  else { 
	    tspk[i][k] = t ; 
	    Volt[i][k] = Vr ; 
	  }
	  
	  if(t>=Tst) { 
	    Mean_Activity[i] += 1. ; 
	    // Idv_Inputs[i][k] += Isyntot[i][k] ; 
	    Idv_Activity[i][k] += 1. ; 
	    if(k<10 && t>=duration+Tst-Tl) fVoltage << " " << Vpeak ; 
	  } 

	  // Updating Postsynaptics inputs 
	  if(IF_INTERPOLATION) {
	    for(j=0;j<nbpop;j++) // i Pres to j Post => Jji 
	      if(J[j][i]!=0) 
	  	for ( l = idxPost[k+Cpt[i]] ; l < (unsigned long) idxPost[k+Cpt[i]] + (unsigned long) nbPost[k+Cpt[i]] ; l++) 
	  	  if(IdPost[l]>=Cpt[j] && IdPost[l]<Cpt[j+1]) 
	  	    Isyn[j][i][IdPost[l]-Cpt[j]] += J[j][i] * exp(-(t-tspk[i][k])/Tsyn[j][i]) ; 
	  }
	  else
	    for(j=0;j<nbpop;j++) // i Pres to j Post => Jji 
	      if(J[j][i]!=0) 
	  	for ( l = idxPost[k+Cpt[i]] ; l < (unsigned long) idxPost[k+Cpt[i]] + (unsigned long) nbPost[k+Cpt[i]] ; l++) 
	  	  if(IdPost[l]>=Cpt[j] && IdPost[l]<Cpt[j+1]) 
	  	    Isyn[j][i][IdPost[l]-Cpt[j]] += J[j][i] ; 
	  
	  //Updating ISI 
	  if(t>=duration+Tst-Tl) 
	    // fprintf(fRaster.c_str(),"%f %f\n",(float) (k+Cpt[i]),(float) (t-Tst) ) ; 
	    fRaster << fixed << setprecision(1) << (float) (k+Cpt[i]) << " " << (float) (tspk[i][k]-duration-Tst+Tl) << endl ; 
	  
	} // endif spike 
	else 
	  if(k<10 && t>=duration+Tst-Tl) fVoltage << " " << Volt[i][k] ;
      } // endfor neurons
    } // endfor populations
    
    if(t>=duration+Tst-Tl) fVoltage << endl ;
    
    // Writing to file  
    // if(t>=Tst && t<=Tst+Tl) {
      
    //   fInput << t-Tst ;
    //   for(int i=0;i<nbpop;i++) { 
    //   	fInput << " " << (Iext[i]+Isyn[i][0][0])*Tm[i] ; 
    //   	for(int j=1;j<nbpop;j++) 
    //   	  fInput << " " << Isyn[i][j][0]*Tm[i] ;
    //   }
    //   fInput << endl ; 
      
    // } //endif 
    
    if(tw>=Tw) {

      cout << " t " << t-Tst-Tw<< " Rates " ;
      for(i=0;i<nbpop;i++) {
    	cout << " " << Mean_Activity[i] / tw*1000. / (double) nbN[i] ; 
	// if(t>Tst+Tl)
	//   Mean_Activity[i] = 0 ; 
      }
      cout.flush();
      
      // if(t<=Tst+Tl) {
	
      // 	fMeanInputs << t-Tst-Tw ; 
      // 	for (int i=0 ;i<nbpop;i++) 
      // 	  for (int j=0 ;j<nbpop;j++) 
      // 	    for (int k=0 ;k<nbN[i];k++) { 
      // 	      fMeanInputs <<" "<< Mean_Isyn[i][j][k] / tw*1000. ; 
      // 	      Mean_Isyn[i][j][k] = 0 ;
      // 	    }	
      // 	fMeanInputs << endl ;
	    
      // 	fIdvInputs << t-Tst-Tw ; 
      // 	for (int i=0 ;i<nbpop;i++) 
      // 	  for (int k=0 ;k<nbN[i];k++) {
      // 	    fIdvInputs <<" "<< Idv_Inputs[i][k] / tw*1000. ;
      // 	    Idv_Inputs[i][k] = 0 ;
      // 	  }
      // 	fIdvInputs << endl ;
      // }

      fMean << t-Tst-Tw ;
      for(i=0;i<nbpop;i++) { 
	fMean << " " << Mean_Activity[i] / tw*1000. / (double) nbN[i] ;
	Mean_Activity[i] = 0 ; 
      }
      fMean << endl ;
      
      fIdvRates << t-Tst-Tw ; 
      for (int i=0 ;i<nbpop;i++) 
      	for (int k=0 ;k<nbN[i];k++) {
      	  fIdvRates <<" "<< Idv_Activity[i][k] / tw*1000. ;
      	  Idv_Activity[i][k] = 0 ;
      	} 
      fIdvRates << endl ; 

      // if(IF_OPSIN)
      // 	External_Input(nbpop,N,nbN,K,Cff,Iext,IextBL,ndI,Jext,L,path) ;

      tw=0;
      
    }//ENDIF
    
    printProgress (percentage) ;
    
    if(IF_TRANSIENT_IEXT)
      tIext += dt ;
    if(tIext>=T_IEXT) { 
      for(i=0;i<nbpop;i++) 
	while(Iext[i]<=IextBL[i]) 
	  Iext[i] += IextBL[i]/nbIext ; 
      tIext = 0 ;
    }

    // updating Time Windows
    if(t>=Tst) {
      tw += dt ;
      
      if(IF_TIMECOURSE) {
	tc += dt ;
	
	if(tc>=Tc) {
	  if(tc>Tc+dTc+Tw) { 
	    if(IF_TRANSIENT) {
	      cout << endl << "Feedforward Off " << endl ; 	      
	      Iext[PrtrPop] -= sqrt(K) * dPrtr * m0 * (Vth-Vr) ;
	      fill(Jext[PrtrPop].begin(), Jext[PrtrPop].end(), 0) ; 
	      IF_TRANSIENT = 0 ; 
	    }
	  }
	  else {
	    if(!IF_TRANSIENT) {	
	      Iext[PrtrPop] += sqrt(K) * dPrtr * m0 * (Vth-Vr) ;
	      cout << endl << "Feedforward On " << PrtrPop << " " << Iext[PrtrPop] << endl ;
	      // if(DIM==1) 
	      // 	External_Input(nbpop,N,nbN,K,Cff,Iext,IextBL,ndI,Jext,path) ;
	      // else
	      // 	External_Input2D(nbpop,N,nbN,K,Cff,Iext,IextBL,ndI,Jext,path) ;
	      IF_TRANSIENT = 1 ;
	    }  
	  }	  
	}
      }
    }
  } //ENDMAINLOOP 

  // fIdvRates << duration ; 
  // for (int i=0 ;i<nbpop;i++) 
  //   for (int k=0 ;k<nbN[i];k++) {
  //     fIdvRates <<" "<< Idv_Activity[i][k] / duration*1000. ;
  //     Idv_Activity[i][k] = 0 ;
  //   } 
  // fIdvRates << endl ; 
  
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
  Mean_Isyn.clear() ; 
  Isyntot.clear() ; 
  Volt.clear() ; 

  Jext.clear() ; 
  Mean_Activity.clear() ; 
  Idv_Inputs.clear() ; 
  Idv_Activity.clear() ; 

  ///////////////////////////////////////////////////////////////////    

  fMean.close();
  // fInput.close();
  fVoltage.close();
  fRaster.close();
  fMeanInputs.close();
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
