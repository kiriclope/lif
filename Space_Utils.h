#ifndef __SPACEUTILS__
#define __SPACEUTILS__

///////////////////////////////////////////////////////////////////////

double DeltaFunc(double X, double Y) {
  if(abs(X-Y)<.00001)
    return 1. ;
  else
    return 0. ;
}

///////////////////////////////////////////////////////////////////    

double Gaussian_1D(double mu, double sigma) {
  return exp(-mu * mu /2.0 /sigma /sigma) / sqrt(2.0*M_PI) /sigma ;
}

///////////////////////////////////////////////////////////////////    

double ShortestDistOnCirc(double point0, double point1) {
  double dist = 0.0;
  dist = abs(point0 - point1) ;
  dist = fmod(dist, L) ;
  if(dist > 0.5*L)
    dist = dist-L ;
  else
    dist = 0.5*L * dist ;
  return dist;
}

///////////////////////////////////////////////////////////////////    

double PeriodicGaussian(double xa, double xb, double varianceOfGaussian) {
  double distX = 0.0; //ShortestDistOnCirc(xa, xb, patchSize);
  distX = ShortestDistOnCirc(xa, xb) ;
  return Gaussian_1D(distX, varianceOfGaussian);
}

///////////////////////////////////////////////////////////////////////

double Wrapped_Gaussian(double mu, double sigma, int klim) {
  double sum = 0 ; 
  for(int k=-klim;k<=klim;k++)
    sum += Gaussian_1D(L*mu+L*(double)k,sigma) ;
  return sum ;
}

///////////////////////////////////////////////////////////////////////

double Wrapped_Exp(double mu, double sigma, int klim) {
  double sum = 0 ; 
  for(int k=-klim;k<=klim;k++)
    sum += exp( -fabs( L* ( mu+(double)k ) ) / sigma ) ;
  return sum ;
}

///////////////////////////////////////////////////////////////////////

void getCrecCff(char** argv, int nbpop, double *&Crec, double &Cff) {

  Crec = new double [nbpop] ; 
  for(int i=0;i<nbpop;i++) 
    Crec[i] = (double) atof(argv[ argIext * nbpop + 6 + i + IF_Prtr ]) ;
  
  if(IF_Prtr)
    Cff = (double) atof(argv[ ( argIext + 1) * nbpop + 6 + IF_Prtr ]) ; // variance of the wrapped external input

}

///////////////////////////////////////////////////////////////////////

void CreateDir_SpaceCrec(int nbpop,string &path,unsigned long N,double* Crec) {
  
  if(PROFILE==1) {
    cout  << "O( sqrt(K) ) Cosine interactions " << endl ;
    if(DIM==1)
      if(IF_SPEC)
	path += "/Ring/Spec/" ;
      else
	path += "/Ring/" ;
    else
      path += "/Ring2D/" ;
  }
  if(PROFILE==2) {
    cout  << "O( sqrt(K) ) Gaussian interactions " << endl ;    
    if(DIM==1)
      path += "/Gauss/" ;
    else
      path += "/Gauss2D/" ;
  }
  if(PROFILE==3) {
    cout  << "O( sqrt(K) ) Exp interactions " << endl ;
    if(DIM==1)
      path += "/Exp/" ;
    else
      path += "/Exp2D/" ;
  }

  string mkdirp = "mkdir -p " ;
  
  string popList[4] = {"E","I","S","V"} ;  
  if(nbpop==1) 
    popList[0] = "I" ;
  string strpop ;
  
  char cCrec[10] ;
  
  for(int i=0;i<nbpop;i++) {
    strpop = popList[i] ;
    sprintf(cCrec,"%0.4f",Crec[i]) ; 
    path += "Crec"+ strpop + string(cCrec) ; 
  }
  mkdirp += path ; 

  const char * cmd = mkdirp.c_str();
  const int dir_err = system(cmd);

  if(-1 == dir_err) {
    cout << "error creating directories" << endl ;
  }
  
  cout << "Created directory : " ;
  cout << path << endl ;
}

///////////////////////////////////////////////////////////////////////

void CreateDir_SpaceCrec2D(int nbpop,string &path,int N,double* Crec) {
  
  if(PROFILE==1) {
    cout  << "O( sqrt(K) ) Cosine interactions " << endl ;
    path += "/Ring2D/" ;
  }
  if(PROFILE==2) {
    cout  << "O( sqrt(K) ) Gaussian interactions " << endl ;
    path += "/Gauss2D/" ;
  }
  string mkdirp = "mkdir -p " ;
  
  string popList[4] = {"E","I","S","V"} ;  
  if(nbpop==1) 
    popList[0] = "I" ;
  string strpop ;
  
  char cCrec[10] ;
  
  for(int i=0;i<nbpop;i++) {
    strpop = popList[i] ;
    sprintf(cCrec,"%0.2f",Crec[i]) ; 
    path += "Crec"+ strpop + string(cCrec) ; 
  }
  mkdirp += path ; 

  const char * cmd = mkdirp.c_str();
  const int dir_err = system(cmd);

  if(-1 == dir_err) {
    cout << "error creating directories" << endl ;
  }
  
  cout << "Created directory : " ;
  cout << path << endl ;
}

///////////////////////////////////////////////////////////////////////

void CreateDir_SpaceCrecCab(int nbpop,string &path,int N,double **Crec) {
  
  if(PROFILE==1) {
    cout  << "O( sqrt(K) ) Cosine interactions " << endl ;
    path += "/Ring/" ;
  }
  if(PROFILE==2) {
    cout  << "O( sqrt(K) ) Gaussian interactions " << endl ;
    path += "/Gauss/" ;
  }
  if(PROFILE==3) {
    cout  << "O( sqrt(K) ) Exp interactions " << endl ;
    path += "/Exp/" ;
  }

  string mkdirp = "mkdir -p " ; 
  string popList[4] = {"E","I","S","V"} ;  
  if(nbpop==1) 
    popList[0] = "I" ;
  string strPre,strPost ;

  char cCrec[10] ;

  for(int i=0;i<nbpop;i++) {
    strPost = popList[i] ;
    // cout << strPost << " " << strPre  << endl ;
    for(int j=0;j<nbpop;j++) {
      strPre = popList[j] ;
      sprintf(cCrec,"%0.4f",Crec[i][j]) ;
      path += "Crec"+ strPost + strPre + string(cCrec) ; 
    }
  }
  mkdirp += path ; 

  const char * cmd = mkdirp.c_str();
  const int dir_err = system(cmd);

  if(-1 == dir_err) {
    cout << "error creating directories" << endl ;
  }
  
  cout << "Created directory : " ;
  cout << path << endl ;
}

void CreateDir_SpaceCff(int nbpop,string &path,int N,double Cff) {
  
  string mkdirp = "mkdir -p " ;
 
  char cCff[10] ;
  sprintf(cCff,"%0.4f",Cff) ;
  string sCff = string(cCff) ;
  
  path += "Cff"+sCff ; 
  mkdirp += path ; 

  const char * cmd = mkdirp.c_str();
  const int dir_err = system(cmd);

  if(-1 == dir_err) {
    cout << "error creating directories" << endl ;
  }
  
  cout << "Created directory : " ;
  cout << path << endl ;
}

///////////////////////////////////////////////////////////////////////

void Spatial_Connection_Probability_1D(int nbpop,int N,int* nbN,vector<int> &Cpt,double K,int klim,double* Crec,double **&c) {

  double **ix ;
  ix = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    ix[i] = new double[nbN[i]] ;

  double **X ;
  X = new double*[nbpop] ; 
  for(int i=0;i<nbpop;i++)
    X[i] = new double[nbN[i]] ;


  //////////////////////////////////////////////////

  for(int i=0;i<nbpop;i++)
    for(int k=0;k<nbN[i];k++) {
      ix[i][k] = fmod( double(k), double( nbN[i]) ) ;
      X[i][k] = ix[i][k]/( double( nbN[i]) ) ;
    }

  cout << "ix " ;
  for(int k=0;k<10;k++)
    cout << ix[0][k] << " " ;
  cout << " " ;
  cout << endl ;
  for(int k=0;k<10;k++)
    cout << ix[0][nbN[0]-1-k] << " " ;
  cout << endl ;

  cout << "X " ;
  for(int k=0;k<10;k++)
    cout << X[0][k] << " " ;
  cout << endl ;
  cout << " " ;
  for(int k=0;k<10;k++)
    cout << X[0][nbN[0]-1-k] << " " ;
  cout << endl ;

  if(nbpop>1) {
    cout << "X " ;
    for(int k=0;k<10;k++)
      cout << X[1][k] << " " ;
    cout << endl ;
    cout << " " ;
    for(int k=0;k<10;k++)
      cout << X[1][nbN[1]-1-k] << " " ;
    cout << endl ;
  }

  double ***sum ;
  sum = new double**[nbpop] ;
  for(int i=0;i<nbpop;i++) {
    sum[i] = new double*[nbpop] ;
    for(int j=0;j<nbpop;j++)
      sum[i][j] = new double[nbN[i]] ;
  }

  c = new double*[N] ; 
  for(int i=0;i<N;i++) 
    c[i] = new double[N] ; 
  
  cout << "Connection Probability : " ;

  if(PROFILE==1) {// O( sqrt(K) ) Cosine interactions
    
    cout  << "O( sqrt(K) ) Cosine interactions " << endl ;
    for(int i=0;i<nbpop;i++)
      for(int k=0;k<nbN[i];k++)
	for(int j=0;j<nbpop;j++) 
	  for(int l=0;l<nbN[j];l++) {
	    c[k+Cpt[i]][l+Cpt[j]] = K / (double) nbN[j] * ( 1.0 + 2.0 * Crec[j] * cos( 2.0 * M_PI * ( X[i][k]-X[j][l] ) ) ) ;
	  }
  }
  
  if(PROFILE==2) { // O( sqrt(K) ) Gaussian interactions
    cout  << "O( sqrt(K) ) Gaussian interactions " << endl ;

    for(int i=0;i<nbpop;i++)
      for(int k=0;k<nbN[i];k++)
    	for(int j=0;j<nbpop;j++) {
    	  for(int l=0;l<nbN[j];l++) {
	    c[k+Cpt[i]][l+Cpt[j]] = Wrapped_Gaussian(X[i][k]-X[j][l],Crec[j],klim) ;
    	    sum[i][j][k] += c[k+Cpt[i]][l+Cpt[j]] ;
	  }
	  for(int l=0;l<nbN[j];l++)
    	    c[k+Cpt[i]][l+Cpt[j]] = K / sum[i][j][k] * c[k+Cpt[i]][l+Cpt[j]] ;
	}
  }

  if(PROFILE==3) { // O( sqrt(K) ) Gaussian interactions
    cout  << "O( sqrt(K) ) Exp interactions " << endl ;

    for(int i=0;i<nbpop;i++)
      for(int k=0;k<nbN[i];k++)
    	for(int j=0;j<nbpop;j++) {
    	  for(int l=0;l<nbN[j];l++) {
	    c[k+Cpt[i]][l+Cpt[j]] = Wrapped_Exp(X[i][k]-X[j][l],Crec[j],klim) ;
    	    sum[i][j][k] += c[k+Cpt[i]][l+Cpt[j]] ;
	  }
	  for(int l=0;l<nbN[j];l++)
    	    c[k+Cpt[i]][l+Cpt[j]] = K / sum[i][j][k] * c[k+Cpt[i]][l+Cpt[j]] ;
	}
  }

  delete [] ix ;
  delete [] X ;
  delete [] sum ;  
}

void Spatial_Connection_Probability_1D_Ka(int nbpop,int N,int* nbN,vector<int> &Cpt,double *K,int klim,double* Crec,double** &c) {

  double **ix ;
  ix = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    ix[i] = new double[nbN[i]] ;

  double **X ;
  X = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    X[i] = new double[nbN[i]] ;

  for(int i=0;i<nbpop;i++)
    for(int k=0;k<nbN[i];k++) {
      ix[i][k] = fmod( double(k), double( nbN[i] ) ) ;
      X[i][k] = ix[i][k]/( double( nbN[i]) ) ;
    }

  cout << "ix " ;
  for(int k=0;k<10;k++) 
    cout << ix[0][k] << " " ;
  cout << endl ;
  for(int k=0;k<10;k++) 
    cout << ix[0][nbN[0]-1-k] << " " ;
  cout << endl ;

  cout << "X " ;
  for(int k=0;k<10;k++) 
    cout << X[0][k] << " " ;
  cout << endl ;
  for(int k=0;k<10;k++) 
    cout << X[0][nbN[0]-1-k] << " " ;
  cout << endl ;

  double ***sum ;
  sum = new double**[nbpop] ;
  for(int i=0;i<nbpop;i++) {
    sum[i] = new double*[nbpop] ;
    for(int j=0;j<nbpop;j++)
      sum[i][j] = new double[nbN[i]] ;
  }

  c = new double*[N] ;      
  for(int i=0;i<N;i++)
    c[i] = new double[N] ;
  
  cout << "Connection Probability ... " << endl ;
  for(int i=0;i<nbpop;i++)
    for(int k=0;k<nbN[i];k++)
      for(int j=0;j<nbpop;j++) 
 	for(int l=0;l<nbN[j];l++) {
	  if(Crec[j]==0)
	    c[k+Cpt[i]][l+Cpt[j]] = K[j]/(double)nbN[j] * DeltaFunc(X[i][k],X[j][l]) ;
	  else {
	    c[k+Cpt[i]][l+Cpt[j]] = Wrapped_Gaussian(X[i][k]-X[j][l],Crec[j],klim) ;
	    sum[i][j][k] = sum[i][j][k] + c[k+Cpt[i]][l+Cpt[j]] ;
	  }
	}

  for(int i=0;i<nbpop;i++)
    for(int k=0;k<nbN[i];k++)
      for(int j=0;j<nbpop;j++)
  	for(int l=0;l<nbN[j];l++)
  	  if(Crec[j]!=0)
  	    c[k+Cpt[i]][l+Cpt[j]] = K[j]/sum[i][j][k]*c[k+Cpt[i]][l+Cpt[j]] ;

  delete [] ix ;
  delete [] X ;
  delete [] sum ;
}

///////////////////////////////////////////////////////////////////////

void Spatial_Connection_Probability_1Dab(int nbpop,int N,int* nbN,vector<int> &Cpt,double K,int klim,double **Crec,double ** &c) {

  double **ix ;
  ix = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    ix[i] = new double[nbN[i]] ;

  double **X ;
  X = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    X[i] = new double[nbN[i]] ;

  for(int i=0;i<nbpop;i++)
    for(int k=0;k<nbN[i];k++) {
      ix[i][k] = fmod( double(k), double( nbN[i]) ) ;
      X[i][k] = ix[i][k]/ (double( nbN[i]) ) ;
    }

  cout << "ix " ;
  for(int k=0;k<10;k++) 
    cout << ix[0][k] << " " ;
  cout << endl ;
  for(int k=0;k<10;k++) 
    cout << ix[0][nbN[0]-1-k] << " " ;
  cout << endl ;

  cout << "X " ;
  for(int k=0;k<10;k++) 
    cout << X[0][k] << " " ;
  cout << endl ;
  for(int k=0;k<10;k++) 
    cout << X[0][nbN[0]-1-k] << " " ;
  cout << endl ;

  double ***sum ;
  sum = new double**[nbpop] ;
  for(int i=0;i<nbpop;i++) {
    sum[i] = new double*[nbpop] ;
    for(int j=0;j<nbpop;j++)
      sum[i][j] = new double[nbN[i]] ;
  }

  c = new double*[N] ;      
  for(int i=0;i<N;i++)
    c[i] = new double[N] ;
  
  cout << "Connection Probability ... " << endl ;
  for(int i=0;i<nbpop;i++)
    for(int k=0;k<nbN[i];k++)
      for(int j=0;j<nbpop;j++) {
 	for(int l=0;l<nbN[j];l++) {
	  c[k+Cpt[i]][l+Cpt[j]] = Wrapped_Gaussian(X[i][k]-X[j][l],Crec[i][j],klim) ;
	  sum[i][j][k] += c[k+Cpt[i]][l+Cpt[j]] ;
	}
	for(int l=0;l<nbN[j];l++) 
	  c[k+Cpt[i]][l+Cpt[j]] = K/sum[i][j][k]*c[k+Cpt[i]][l+Cpt[j]] ;
      }
  
  delete [] ix ;
  delete [] X ;
  delete [] sum ;  
}

//////////////////////////////////////////////////////////////////


void Spatial_Connection_Probability_2D(int nbpop,int N,int* nbN,vector<int> &Cpt,double K,int klim,double* Crec,double** &c) {

  double **ix ;
  ix = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    ix[i] = new double[nbN[i]] ;

  double **iy ;
  iy = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    iy[i] = new double[nbN[i]] ;

  double **X ;
  X = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    X[i] = new double[nbN[i]] ;

  double **Y ;
  Y = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    Y[i] = new double[nbN[i]] ;

  for(int i=0;i<nbpop;i++)
    for(int k=0;k<nbN[i];k++) {
      ix[i][k] = fmod( double(k), sqrt( double( nbN[i]) ) ) ;
      X[i][k] = ix[i][k]/sqrt( double( nbN[i]) ) ;

      iy[i][k] = floor( double(k)/sqrt( double( nbN[i]) ) ) ;
      Y[i][k] = iy[i][k]/sqrt( double (nbN[i]) ) ;
    }

  cout << "ix " ;
  for(int k=0;k<10;k++) 
    cout << ix[0][k] << " " ;
  cout << endl ;
  for(int k=0;k<10;k++) 
    cout << ix[0][nbN[0]-1-k] << " " ;
  cout << endl ;

  cout << "iy " ;
  for(int k=0;k<10;k++) 
    cout << iy[0][k] << " " ;
  cout << endl ;
  for(int k=0;k<10;k++) 
    cout << iy[0][nbN[0]-1-k] << " " ;
  cout << endl ;

  cout << "X " ;
  for(int k=0;k<10;k++) 
    cout << X[0][k] << " " ;
  cout << endl ;
  for(int k=0;k<10;k++) 
    cout << X[0][nbN[0]-1-k] << " " ;
  cout << endl ;

  cout << "Y " ;
  for(int k=0;k<10;k++) 
    cout << Y[0][k] << " " ;
  cout << endl ;
  for(int k=0;k<10;k++) 
    cout << Y[0][nbN[0]-1-k] << " " ;
  cout << endl ;
  
  ///////////////////////////////////////////////////////////////////////
  // Spatiality 
  ///////////////////////////////////////////////////////////////////////

  double ***sum ;
  sum = new double**[nbpop] ;
  for(int i=0;i<nbpop;i++) {
    sum[i] = new double*[nbpop] ;
    for(int j=0;j<nbpop;j++)
      sum[i][j] = new double[nbN[i]] ;
  }

  c = new double*[N] ;      
  for(int i=0;i<N;i++)
    c[i] = new double[N] ;
  
  cout << "Connection Probability ... " << endl ;

  if(PROFILE==1) {// O(1) Cosine interactions
    cout  << "O(1) toroidal interactions " << endl ;
    for(int i=0;i<nbpop;i++)
      for(int k=0;k<nbN[i];k++)
	for(int j=0;j<nbpop;j++) 
	  for(int l=0;l<nbN[j];l++)
	    if(i==0 & j==0)
	      c[k+Cpt[i]][l+Cpt[j]] = K / (double) nbN[j] * ( 1.0 + 2.0 * Crec[j] / sqrt(K) * ( cos( 2 * M_PI * ( X[i][k]-X[j][l] ) ) +  cos( 2 * M_PI * ( Y[i][k]-Y[j][l] ) ) ) ) ;
	    else					 
	      c[k+Cpt[i]][l+Cpt[j]] = K / (double) nbN[j] ;
  }
  if(PROFILE==2) { // O(1) Gaussian interactions
    cout  << "O(1) Gaussian interactions " << endl ;

    for(int i=0;i<1;i++)
      for(int k=0;k<nbN[i];k++)
	for(int j=0;j<1;j++) 
	  for(int l=0;l<nbN[j];l++) {
	    c[k+Cpt[i]][l+Cpt[j]] = Wrapped_Gaussian(X[i][k]-X[j][l],Crec[j],klim)*Wrapped_Gaussian(Y[i][k]-Y[j][l],Crec[j],klim) ;
	    sum[i][j][k] += c[k+Cpt[i]][l+Cpt[j]] ;
	  }

    for(int i=0;i<nbpop;i++)
      for(int k=0;k<nbN[i];k++)
	for(int j=0;j<nbpop;j++)
	  for(int l=0;l<nbN[j];l++)
	    if(i==0 & j==0)
	      c[k+Cpt[i]][l+Cpt[j]] = K / (double) nbN[j] + sqrt(K) / sum[i][j][k] * c[k+Cpt[i]][l+Cpt[j]] ;
	    else
	      c[k+Cpt[i]][l+Cpt[j]] = K / (double) nbN[j] ;
  }
  
  delete [] ix ;
  delete [] iy ;
  delete [] X ;
  delete [] Y ;
  delete [] sum ;
}

///////////////////////////////////////////////////////////////////////

void External_Input(int nbpop, unsigned long N, unsigned long* nbN, double K, double Cff, double* Iext, double *IextBL, int ndI, vector<vector<double> > &Jext, string path) {

  cout << "External Input : " ;

  double **X ;
  X = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    X[i] = new double[nbN[i]] ;

  for(int i=0;i<nbpop;i++)
    for(unsigned long j=0;j<nbN[i];j++) 
      X[i][j] = L * fmod( double(j), double( nbN[i]) - 1.0) / ( double( nbN[i] ) - 1.0 ) ;
    
  double p = 1.0 ;
  
  cout << "Baseline " ;
  for(int i=0;i<nbpop;i++)
    cout << IextBL[i]/sqrt(K)/m0/(Vth-Vr) << " " ;

  if(abs(Iext[ndI]-IextBL[ndI])<=.001) {
    if(IF_GAUSS==-1) {
      cout << "| FeedForward On " << endl ;       
      for(int i=0;i<nbpop;i++)
	for(unsigned long j=0;j<nbN[i];j++) 
	  Jext[i][j] = IextBL[i] * Cff * cos( 2. * ( X[i][j] - X[i][nbN[i]/2-1] - PHI0 ) ) ; 
    }
    else
      cout << endl ;
  }
  else {
    p = ( Iext[ndI] - IextBL[ndI] ) ; 
    
    cout << "| Perturbation " ;
    for(int i=0;i<nbpop;i++)
      cout << Iext[i]/sqrt(K)/m0/(Vth-Vr) - IextBL[i]/sqrt(K)/m0/(Vth-Vr) << " " ;
    cout << endl ;

    for(unsigned long j=0;j<nbN[ndI];j++) { 
      
      switch(IF_GAUSS) {
	
      case 0 :
	if( abs( X[ndI][j]-X[ndI][nbN[ndI]/2-1] ) <= Cff / 2.0 ) 
	  Jext[ndI][j] = p ; 
	else 
	  Jext[ndI][j] = 0. ;  
	break ; 
	
      case 1 :
	
	if( abs( X[ndI][j]-X[ndI][nbN[ndI]/2-1] ) <= 4.0 * Cff )
	  Jext[ndI][j] = p*PeriodicGaussian(X[ndI][j],X[ndI][nbN[ndI]/2-1],Cff) / sqrt(K) ;
	else
	  Jext[ndI][j] = 0. ;
	break ;

      case 2 :

	if( abs(X[ndI][j]-X[ndI][nbN[ndI]/2-1] ) <= Cff )
	  Jext[ndI][j] = p ; // half height of the Gaussian
	else
	  if( abs(X[ndI][j]-X[ndI][nbN[ndI]/2-1] ) <= 4.0 * Cff ) 
	    Jext[ndI][j] = p*PeriodicGaussian(X[ndI][j],X[ndI][nbN[ndI]/2-1],Cff) / ( exp(- 1.0 / 2.0) / sqrt(2.0*M_PI) / Cff ) ; 
	  else 
	    Jext[ndI][j] = 0. ; 
	break ;

      default :
	if( abs( X[ndI][j]-X[ndI][nbN[ndI]/2-1] ) <= Cff / 2.0 ) 
	  Jext[ndI][j] = p ; 
	else
	  Jext[ndI][j] = 0. ; 

	break ; 

      }	
      
    }
    
  }
    
  /* random_device rd ; */
  /* // default_random_engine gen( rd() ) ; */
  /* default_random_engine gen( 123456789 ) ; */
  /* uniform_real_distribution<double> unif(0,1) ; */

  /* /\* cout << "Check random seed " ; *\/ */
  /* /\* for(int i=0;i<10;i++) *\/ */
  /* /\*   cout << unif(gen) << " " ; *\/ */
  /* /\* cout << endl ; *\/ */

  /* double Pb=1. ; */
  /* if(IF_OPSIN) */
  /*   Pb = OpsPb ; */

  /* string strIdxFile = path + "/IdxFile.txt" ; */
  /* ofstream fIdxFile(strIdxFile.c_str(), ios::out | ios::ate); */
    
  delete [] X ;
}

///////////////////////////////////////////////////////////////////////

void External_Input2D(int nbpop, unsigned long N, unsigned long* nbN, double K, double Cff, double* Iext, double *IextBL, int ndI, vector<vector<double> > &Jext, string path) {

  cout << "External Input : " ;

  double **X ;
  X = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    X[i] = new double[nbN[i]] ;

  for(int i=0;i<nbpop;i++)
    for(unsigned long j=0;j<nbN[i];j++) 
      X[i][j] = L * fmod( double(j), sqrt( double( nbN[i]) ) ) / sqrt( double( nbN[i] ) ) ;
  
  double **Y ;
  Y = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    Y[i] = new double[nbN[i]] ;

  for(int i=0;i<nbpop;i++)
    for(unsigned long j=0;j<nbN[i];j++)
      Y[i][j] = L * floor( double(j) / sqrt( double( nbN[i]) ) ) / sqrt( double (nbN[i]) ) ;
  
  double p = 1.0 ;
  
  cout << "Baseline " ;
  for(int i=0;i<nbpop;i++)
    cout << IextBL[i]/sqrt(K)/m0 << " " ;
  
  if(abs(Iext[ndI]-IextBL[ndI])<=.001) {
    cout << endl ;
  }
  else {
    p = Iext[ndI] - IextBL[ndI] ;
    
    cout << "| Perturbation " ;
    for(int i=0;i<nbpop;i++)
      cout << Iext[i]/sqrt(K)/m0 - IextBL[i]/sqrt(K)/m0 << " " ;
    cout << endl ;

    for(unsigned long j=0;j<nbN[ndI];j++) {

      switch(IF_GAUSS) {
	
      case 0 :
	if( ( X[ndI][j]-.5*L ) * ( X[ndI][j]-.5*L ) + ( Y[ndI][j]-.5*L ) * ( Y[ndI][j]-.5*L ) <= Cff*Cff / 4.0 ) 
	  Jext[ndI][j] = p ; 
	else
	  Jext[ndI][j] = 0. ; 
	break ;
      
      case 1 :
	
	if( abs( X[ndI][j]-X[ndI][nbN[ndI]/2-1] ) <= 4.0 * Cff )
	    Jext[ndI][j] = p*PeriodicGaussian(X[ndI][j],X[ndI][nbN[ndI]/2-1],Cff) ;
	  else
	    Jext[ndI][j] = 0. ;
	break ;

      case 2 :

	if( abs(X[ndI][j]-X[ndI][nbN[ndI]/2-1] ) <= Cff )
	  Jext[ndI][j] = p ; // half height of the Gaussian
	else
	  if( abs(X[ndI][j]-X[ndI][nbN[ndI]/2-1] ) <= 4.0 * Cff )
	    Jext[ndI][j] = p*PeriodicGaussian(X[ndI][j],X[ndI][nbN[ndI]/2-1],Cff) / ( exp(- 1.0 / 2.0) / sqrt(2.0*M_PI) / Cff ) ;
	  else
	    Jext[ndI][j] = 0. ;
	break ;

      default :
	if( abs( X[ndI][j]-X[ndI][nbN[ndI]/2-1] ) <= Cff / 2.0 ) 
	  Jext[ndI][j] = p ; 
	else
	  Jext[ndI][j] = 0. ; 

	break ;

      }	
      
    }
    
  }
        
  delete [] X ;
  delete [] Y ;
}

//////////////////////////////////////////////////////////////////

#endif
