#ifndef __SPACEUTILS__
#define __SPACEUTILS__

#define PROFILE 2
#define BL .25
#define m0 1E-2

#define IF_OPSIN false
#define OpsPb 1.

///////////////////////////////////////////////////////////////////////

double DeltaFunc(double X, double Y) {
  if(abs(X-Y)<.00001)
    return 1. ;
  else
    return 0. ;
}

double Gaussian_1D(double mu, double sigma) {
  return exp(-mu*mu/2./sigma/sigma)/sqrt(2.*M_PI)/sigma ;
}

///////////////////////////////////////////////////////////////////////

double Wrapped_Gaussian(double mu, double sigma, int klim, double L) {
  double sum = 0 ; 
  for(int k=-klim;k<=klim;k++)
    sum += Gaussian_1D(L*mu+L*(double)k,sigma) ;
  return sum ;
}

///////////////////////////////////////////////////////////////////////

double Wrapped_Exp(double mu, double sigma, int klim, double L) {
  double sum = 0 ; 
  for(int k=-klim;k<=klim;k++)
    sum += exp( -fabs( L* ( mu+(double)k ) ) / sigma ) ;
  return sum ;
}

///////////////////////////////////////////////////////////////////////

void CreateDir_SpaceCrec(int nbpop,string &path,int N,double* Crec) {
  
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

void Spatial_Connection_Probability_1D(int nbpop,int N,int* Nk,vector<int> &Cpt,double K,int klim,double* Crec,double **&c,double L) {

  double **ix ;
  ix = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    ix[i] = new double[Nk[i]] ;

  double **X ;
  X = new double*[nbpop] ; 
  for(int i=0;i<nbpop;i++)
    X[i] = new double[Nk[i]] ;


  //////////////////////////////////////////////////

  for(int i=0;i<nbpop;i++)
    for(int k=0;k<Nk[i];k++) {
      ix[i][k] = fmod( double(k), double( Nk[i]) ) ;
      X[i][k] = ix[i][k]/( double( Nk[i]) ) ;
    }

  cout << "ix " ;
  for(int k=0;k<10;k++)
    cout << ix[0][k] << " " ;
  cout << " " ;
  cout << endl ;
  for(int k=0;k<10;k++)
    cout << ix[0][Nk[0]-1-k] << " " ;
  cout << endl ;

  cout << "X " ;
  for(int k=0;k<10;k++)
    cout << X[0][k] << " " ;
  cout << endl ;
  cout << " " ;
  for(int k=0;k<10;k++)
    cout << X[0][Nk[0]-1-k] << " " ;
  cout << endl ;

  if(nbpop>1) {
    cout << "X " ;
    for(int k=0;k<10;k++)
      cout << X[1][k] << " " ;
    cout << endl ;
    cout << " " ;
    for(int k=0;k<10;k++)
      cout << X[1][Nk[1]-1-k] << " " ;
    cout << endl ;
  }

  double ***sum ;
  sum = new double**[nbpop] ;
  for(int i=0;i<nbpop;i++) {
    sum[i] = new double*[nbpop] ;
    for(int j=0;j<nbpop;j++)
      sum[i][j] = new double[Nk[i]] ;
  }

  c = new double*[N] ; 
  for(int i=0;i<N;i++) 
    c[i] = new double[N] ; 
  
  cout << "Connection Probability : " ;

  if(PROFILE==1) {// O( sqrt(K) ) Cosine interactions
    
    cout  << "O( sqrt(K) ) Cosine interactions " << endl ;
    for(int i=0;i<nbpop;i++)
      for(int k=0;k<Nk[i];k++)
	for(int j=0;j<nbpop;j++) 
	  for(int l=0;l<Nk[j];l++) {
	    c[k+Cpt[i]][l+Cpt[j]] = K / (double) Nk[j] * ( 1.0 + 2.0 * Crec[j] * cos( 2.0 * M_PI * ( X[i][k]-X[j][l] ) ) ) ;
	  }
  }
  
  if(PROFILE==2) { // O( sqrt(K) ) Gaussian interactions
    cout  << "O( sqrt(K) ) Gaussian interactions " << endl ;

    for(int i=0;i<nbpop;i++)
      for(int k=0;k<Nk[i];k++)
    	for(int j=0;j<nbpop;j++) {
    	  for(int l=0;l<Nk[j];l++) {
	    c[k+Cpt[i]][l+Cpt[j]] = Wrapped_Gaussian(X[i][k]-X[j][l],Crec[j],klim,L) ;
    	    sum[i][j][k] += c[k+Cpt[i]][l+Cpt[j]] ;
	  }
	  for(int l=0;l<Nk[j];l++)
    	    c[k+Cpt[i]][l+Cpt[j]] = K / sum[i][j][k] * c[k+Cpt[i]][l+Cpt[j]] ;
	}
  }

  if(PROFILE==3) { // O( sqrt(K) ) Gaussian interactions
    cout  << "O( sqrt(K) ) Exp interactions " << endl ;

    for(int i=0;i<nbpop;i++)
      for(int k=0;k<Nk[i];k++)
    	for(int j=0;j<nbpop;j++) {
    	  for(int l=0;l<Nk[j];l++) {
	    c[k+Cpt[i]][l+Cpt[j]] = Wrapped_Exp(X[i][k]-X[j][l],Crec[j],klim,L) ;
    	    sum[i][j][k] += c[k+Cpt[i]][l+Cpt[j]] ;
	  }
	  for(int l=0;l<Nk[j];l++)
    	    c[k+Cpt[i]][l+Cpt[j]] = K / sum[i][j][k] * c[k+Cpt[i]][l+Cpt[j]] ;
	}
  }

  delete [] ix ;
  delete [] X ;
  delete [] sum ;  
}

void Spatial_Connection_Probability_1D_Ka(int nbpop,int N,int* Nk,vector<int> &Cpt,double *K,int klim,double* Crec,double** &c,double L) {

  double **ix ;
  ix = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    ix[i] = new double[Nk[i]] ;

  double **X ;
  X = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    X[i] = new double[Nk[i]] ;

  for(int i=0;i<nbpop;i++)
    for(int k=0;k<Nk[i];k++) {
      ix[i][k] = fmod( double(k), double( Nk[i] ) ) ;
      X[i][k] = ix[i][k]/( double( Nk[i]) ) ;
    }

  cout << "ix " ;
  for(int k=0;k<10;k++) 
    cout << ix[0][k] << " " ;
  cout << endl ;
  for(int k=0;k<10;k++) 
    cout << ix[0][Nk[0]-1-k] << " " ;
  cout << endl ;

  cout << "X " ;
  for(int k=0;k<10;k++) 
    cout << X[0][k] << " " ;
  cout << endl ;
  for(int k=0;k<10;k++) 
    cout << X[0][Nk[0]-1-k] << " " ;
  cout << endl ;

  double ***sum ;
  sum = new double**[nbpop] ;
  for(int i=0;i<nbpop;i++) {
    sum[i] = new double*[nbpop] ;
    for(int j=0;j<nbpop;j++)
      sum[i][j] = new double[Nk[i]] ;
  }

  c = new double*[N] ;      
  for(int i=0;i<N;i++)
    c[i] = new double[N] ;
  
  cout << "Connection Probability ... " << endl ;
  for(int i=0;i<nbpop;i++)
    for(int k=0;k<Nk[i];k++)
      for(int j=0;j<nbpop;j++) 
 	for(int l=0;l<Nk[j];l++) {
	  if(Crec[j]==0)
	    c[k+Cpt[i]][l+Cpt[j]] = K[j]/(double)Nk[j] * DeltaFunc(X[i][k],X[j][l]) ;
	  else {
	    c[k+Cpt[i]][l+Cpt[j]] = Wrapped_Gaussian(X[i][k]-X[j][l],Crec[j],klim,L) ;
	    sum[i][j][k] = sum[i][j][k] + c[k+Cpt[i]][l+Cpt[j]] ;
	  }
	}

  for(int i=0;i<nbpop;i++)
    for(int k=0;k<Nk[i];k++)
      for(int j=0;j<nbpop;j++)
  	for(int l=0;l<Nk[j];l++)
  	  if(Crec[j]!=0)
  	    c[k+Cpt[i]][l+Cpt[j]] = K[j]/sum[i][j][k]*c[k+Cpt[i]][l+Cpt[j]] ;

  delete [] ix ;
  delete [] X ;
  delete [] sum ;
}

///////////////////////////////////////////////////////////////////////

void Spatial_Connection_Probability_1Dab(int nbpop,int N,int* Nk,vector<int> &Cpt,double K,int klim,double **Crec,double ** &c,double L) {

  double **ix ;
  ix = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    ix[i] = new double[Nk[i]] ;

  double **X ;
  X = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    X[i] = new double[Nk[i]] ;

  for(int i=0;i<nbpop;i++)
    for(int k=0;k<Nk[i];k++) {
      ix[i][k] = fmod( double(k), double( Nk[i]) ) ;
      X[i][k] = ix[i][k]/ (double( Nk[i]) ) ;
    }

  cout << "ix " ;
  for(int k=0;k<10;k++) 
    cout << ix[0][k] << " " ;
  cout << endl ;
  for(int k=0;k<10;k++) 
    cout << ix[0][Nk[0]-1-k] << " " ;
  cout << endl ;

  cout << "X " ;
  for(int k=0;k<10;k++) 
    cout << X[0][k] << " " ;
  cout << endl ;
  for(int k=0;k<10;k++) 
    cout << X[0][Nk[0]-1-k] << " " ;
  cout << endl ;

  double ***sum ;
  sum = new double**[nbpop] ;
  for(int i=0;i<nbpop;i++) {
    sum[i] = new double*[nbpop] ;
    for(int j=0;j<nbpop;j++)
      sum[i][j] = new double[Nk[i]] ;
  }

  c = new double*[N] ;      
  for(int i=0;i<N;i++)
    c[i] = new double[N] ;
  
  cout << "Connection Probability ... " << endl ;
  for(int i=0;i<nbpop;i++)
    for(int k=0;k<Nk[i];k++)
      for(int j=0;j<nbpop;j++) {
 	for(int l=0;l<Nk[j];l++) {
	  c[k+Cpt[i]][l+Cpt[j]] = Wrapped_Gaussian(X[i][k]-X[j][l],Crec[i][j],klim,L) ;
	  sum[i][j][k] += c[k+Cpt[i]][l+Cpt[j]] ;
	}
	for(int l=0;l<Nk[j];l++) 
	  c[k+Cpt[i]][l+Cpt[j]] = K/sum[i][j][k]*c[k+Cpt[i]][l+Cpt[j]] ;
      }
  
  delete [] ix ;
  delete [] X ;
  delete [] sum ;  
}

//////////////////////////////////////////////////////////////////


void Spatial_Connection_Probability_2D(int nbpop,int N,int* Nk,vector<int> &Cpt,double K,int klim,double* Crec,double** &c,double L) {

  double **ix ;
  ix = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    ix[i] = new double[Nk[i]] ;

  double **iy ;
  iy = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    iy[i] = new double[Nk[i]] ;

  double **X ;
  X = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    X[i] = new double[Nk[i]] ;

  double **Y ;
  Y = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    Y[i] = new double[Nk[i]] ;

  for(int i=0;i<nbpop;i++)
    for(int k=0;k<Nk[i];k++) {
      ix[i][k] = fmod( double(k), sqrt( double( Nk[i]) ) ) ;
      X[i][k] = ix[i][k]/sqrt( double( Nk[i]) ) ;

      iy[i][k] = floor( double(k)/sqrt( double( Nk[i]) ) ) ;
      Y[i][k] = iy[i][k]/sqrt( double (Nk[i]) ) ;
    }

  cout << "ix " ;
  for(int k=0;k<10;k++) 
    cout << ix[0][k] << " " ;
  cout << endl ;
  for(int k=0;k<10;k++) 
    cout << ix[0][Nk[0]-1-k] << " " ;
  cout << endl ;

  cout << "iy " ;
  for(int k=0;k<10;k++) 
    cout << iy[0][k] << " " ;
  cout << endl ;
  for(int k=0;k<10;k++) 
    cout << iy[0][Nk[0]-1-k] << " " ;
  cout << endl ;

  cout << "X " ;
  for(int k=0;k<10;k++) 
    cout << X[0][k] << " " ;
  cout << endl ;
  for(int k=0;k<10;k++) 
    cout << X[0][Nk[0]-1-k] << " " ;
  cout << endl ;

  cout << "Y " ;
  for(int k=0;k<10;k++) 
    cout << Y[0][k] << " " ;
  cout << endl ;
  for(int k=0;k<10;k++) 
    cout << Y[0][Nk[0]-1-k] << " " ;
  cout << endl ;
  
  ///////////////////////////////////////////////////////////////////////
  // Spatiality 
  ///////////////////////////////////////////////////////////////////////

  double ***sum ;
  sum = new double**[nbpop] ;
  for(int i=0;i<nbpop;i++) {
    sum[i] = new double*[nbpop] ;
    for(int j=0;j<nbpop;j++)
      sum[i][j] = new double[Nk[i]] ;
  }

  c = new double*[N] ;      
  for(int i=0;i<N;i++)
    c[i] = new double[N] ;
  
  cout << "Connection Probability ... " << endl ;

  if(PROFILE==1) {// O(1) Cosine interactions
    cout  << "O(1) toroidal interactions " << endl ;
    for(int i=0;i<nbpop;i++)
      for(int k=0;k<Nk[i];k++)
	for(int j=0;j<nbpop;j++) 
	  for(int l=0;l<Nk[j];l++)
	    if(i==0 & j==0)
	      c[k+Cpt[i]][l+Cpt[j]] = K / (double) Nk[j] * ( 1.0 + 2.0 * Crec[j] / sqrt(K) * ( cos( 2 * M_PI * ( X[i][k]-X[j][l] ) ) +  cos( 2 * M_PI * ( Y[i][k]-Y[j][l] ) ) ) ) ;
	    else					 
	      c[k+Cpt[i]][l+Cpt[j]] = K / (double) Nk[j] ;
  }
  if(PROFILE==2) { // O(1) Gaussian interactions
    cout  << "O(1) Gaussian interactions " << endl ;

    for(int i=0;i<1;i++)
      for(int k=0;k<Nk[i];k++)
	for(int j=0;j<1;j++) 
	  for(int l=0;l<Nk[j];l++) {
	    c[k+Cpt[i]][l+Cpt[j]] = Wrapped_Gaussian(X[i][k]-X[j][l],Crec[j],klim,L)*Wrapped_Gaussian(Y[i][k]-Y[j][l],Crec[j],klim,L) ;
	    sum[i][j][k] += c[k+Cpt[i]][l+Cpt[j]] ;
	  }

    for(int i=0;i<nbpop;i++)
      for(int k=0;k<Nk[i];k++)
	for(int j=0;j<nbpop;j++)
	  for(int l=0;l<Nk[j];l++)
	    if(i==0 & j==0)
	      c[k+Cpt[i]][l+Cpt[j]] = K / (double) Nk[j] + sqrt(K) / sum[i][j][k] * c[k+Cpt[i]][l+Cpt[j]] ;
	    else
	      c[k+Cpt[i]][l+Cpt[j]] = K / (double) Nk[j] ;
  }
  
  delete [] ix ;
  delete [] iy ;
  delete [] X ;
  delete [] Y ;
  delete [] sum ;
}

///////////////////////////////////////////////////////////////////////

void External_Input(int nbpop,int N,int* Nk,double K,double Cff,double* Iext,double *IextBL,int ndI,vector<vector<double> > &Jext,double L,string path) {

  cout << "External Input" << endl ;

  double **ix ;
  ix = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    ix[i] = new double[Nk[i]] ;

  double **X ;
  X = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    X[i] = new double[Nk[i]] ;

  for(int i=0;i<nbpop;i++)
    for(int k=0;k<Nk[i];k++) {
      ix[i][k] = fmod( double(k), double( Nk[i]) ) ;
      X[i][k] = ix[i][k]/( double( Nk[i] ) ) ;
    }

  cout << "ix " ;
  for(int k=0;k<10;k++) 
    cout << ix[0][k] << " " ;
  cout << endl ;
  for(int k=0;k<10;k++) 
    cout << ix[0][Nk[0]-1-k] << " " ;
  cout << endl ;

  cout << "X " ;
  for(int k=0;k<10;k++) 
    cout << X[0][k] << " " ;
  cout << endl ;
  for(int k=0;k<10;k++) 
    cout << X[0][Nk[0]-1-k] << " " ;
  cout << endl ;
  
  double p=1. ;
  bool BASELINE=0 ;
  if(abs(Iext[ndI]-IextBL[ndI])<=.001) {
    cout << "Baseline => No Input Modulation " << endl;
    BASELINE=1 ;
  }
  else
    cout << "Perturbation" << endl ;

  p = Iext[ndI] - IextBL[ndI] ;

  cout << "Baseline " ;
  for(int i=0;i<nbpop;i++)
    cout << IextBL[i]/sqrt(K)/m0 << " " ;
  cout << endl ;

  if(BASELINE==0)  {
    cout << "Perturbation " ;
    for(int i=0;i<nbpop;i++)
      cout << Iext[i]/sqrt(K)/m0 - IextBL[i]/sqrt(K)/m0 << " " ;
  }

  random_device rd ;
  // default_random_engine gen( rd() ) ;
  default_random_engine gen( 123456789 ) ;
  uniform_real_distribution<double> unif(0,1) ;

  cout << "Check random seed " ;
  for(int i=0;i<10;i++)
    cout << unif(gen) << " " ;
  cout << endl ;

  double Pb=1. ;
  if(IF_OPSIN)
    Pb = OpsPb ;


  string strIdxFile = path + "/IdxFile.txt" ;
  ofstream fIdxFile(strIdxFile.c_str(), ios::out | ios::ate);

  for(int i=0;i<nbpop;i++) {
    for(int j=0;j<Nk[i];j++) {
      
      /* if(i==ndI && BASELINE==0) // Modulated input only onto I */
      /* Jext[i][j] = p*Wrapped_Gaussian(X[i][j]-X[i][Nk[i]/2-1],Cff,100,L) ; */
      /* if(i==ndI && BASELINE==0 && abs(2.*M_PI*(X[i][j]-X[i][Nk[i]/2]))<=Cff) // RectIn */
      /* 	Jext[i][j] = p ; */
      
      /* if(i==ndI && BASELINE==0) // Modulated input only onto I */
      /* if( i==ndI && BASELINE==0 && abs(X[i][j]-X[i][Nk[i]/2-1])<=Cff/L && unif(gen)<=5 ) // Truncated Gaussian */
      if( i==ndI && BASELINE==0 && abs(X[i][j]-X[i][Nk[i]/2-1])<=Cff/L/2 && unif(gen)<=Pb ) {// Truncated Gaussian
	if(IF_OPSIN)
	  Jext[i][j] = p ;
	else
	  Jext[i][j] = p ;
    	fIdxFile << j << " " ;
	/* Jext[i][j] = sqrt(K)*p*Wrapped_Gaussian(X[i][j]-X[i][Nk[i]/2-1],Cff,100,L) ; */
	/* Jext[i][j] = p*Wrapped_Gaussian(X[i][j]-X[i][Nk[i]/2],Cff,100,L) ; */
      } 
      else
	Jext[i][j] = 0. ;
    }
  }
  
  fIdxFile.close() ;
  
  if( BASELINE==0 ) {
    cout << "Jext " ;
    for(int i=0;i<nbpop;i++) {
      for(int j=0;j<10;j++) 
	cout << Jext[i][j] << " " ;
      cout << endl ;
    }
  }

  delete [] ix ;
  delete [] X ;
}

//////////////////////////////////////////////////////////////////

void External_Input_Ka(int nbpop,int N,int* Nk,double *K,double Cff,double* Iext,int ndI,vector<vector<double> > &Jext,double L) {

  cout << "External Input" << endl ;

  double **ix ;
  ix = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    ix[i] = new double[Nk[i]] ;

  double **X ;
  X = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    X[i] = new double[Nk[i]] ;

  for(int i=0;i<nbpop;i++)
    for(int k=0;k<Nk[i];k++) {
      ix[i][k] = fmod( double(k), double( Nk[i]) ) ;
      X[i][k] = ix[i][k]/( double( Nk[i]) ) ;
    }

  cout << "ix " ;
  for(int k=0;k<10;k++) 
    cout << ix[0][k] << " " ;
  cout << endl ;
  for(int k=0;k<10;k++) 
    cout << ix[0][Nk[0]-1-k] << " " ;
  cout << endl ;

  cout << "X " ;
  for(int k=0;k<10;k++) 
    cout << X[0][k] << " " ;
  cout << endl ;
  for(int k=0;k<10;k++) 
    cout << X[0][Nk[0]-1-k] << " " ;
  cout << endl ;
  
  double p=1. ;
  bool BASELINE=0 ;
  if(abs(Iext[ndI]-m0*sqrt(K[ndI])*BL)<=m0) {
    cout << "Baseline => No Input Modulation : " << endl ;
    BASELINE=1 ;
  }

  vector<double> IextBL(nbpop) ;
  if(nbpop==1)
    IextBL[0] = BL*sqrt(K[0])*m0 ;
  else{
    IextBL[0] = Iext[0] ;
    IextBL[1] = BL*sqrt(K[1])*m0 ;
  }

  p = Iext[ndI] - IextBL[ndI] ;

  cout << "Baseline " ;
  for(int i=0;i<nbpop;i++)
    cout << IextBL[i]/sqrt(K[i])/m0 << " " ;
  cout << endl ;

  if(BASELINE==0)  {
    cout << "Perturbation " ;
    for(int i=0;i<nbpop;i++)
      cout << Iext[i]/sqrt(K[i])/m0 << " " ;
  }

  random_device rd ;
  default_random_engine gen( rd() ) ;
  uniform_real_distribution<double> unif(0,1) ;

  for(int i=0;i<nbpop;i++) {
    for(int j=0;j<Nk[i];j++) {      
      if(i==ndI && BASELINE==0) // Modulated input only onto I
      	Jext[i][j] = p*Wrapped_Gaussian(X[i][j]-X[i][Nk[i]/2-1],Cff,100,L) ;
      /* if(i==ndI && BASELINE==0 && abs(X[i][j]-X[i][Nk[i]/2-1])<=Cff/L) // RectIn */
      /* 	Jext[i][j] = p ; */
      /* if(i==ndI && BASELINE==0 && abs(X[i][j]-X[i][Nk[i]/2-1])<=4.*Cff/L && unif(gen)<=5 ) // Truncated Gaussian */
      /* 	Jext[i][j] = p*Wrapped_Gaussian(X[i][j]-X[i][Nk[i]/2-1],Cff,100,L) ; */
      else
	Jext[i][j] = 0. ;
    }
  }
  
  for(int i=0;i<nbpop;i++) {
    for(int j=0;j<10;j++) 
      cout << Jext[i][j] << " " ;
    cout << endl ;
  }
  
  delete [] ix ;
  delete [] X ;
}

///////////////////////////////////////////////////////////////////////

void External_Input2D(int nbpop,int N,int* Nk,double K,double Cff,double* Iext,int ndI,vector<vector<double> > &Jext,double L) {

  cout << "External Input" << endl ;

  double **ix ;
  ix = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    ix[i] = new double[Nk[i]] ;

  double **iy ;
  iy = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    iy[i] = new double[Nk[i]] ;

  double **X ;
  X = new double*[nbpop] ;
  for(int i=0;i<nbpop;i++)
    X[i] = new double[Nk[i]] ;

  double **Y ;
  Y = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    Y[i] = new double[Nk[i]] ;

  for(int i=0;i<nbpop;i++)
    for(int k=0;k<Nk[i];k++) {
      ix[i][k] = fmod( double(k), sqrt( double( Nk[i]) ) ) ;
      X[i][k] = ix[i][k]/sqrt( double( Nk[i]) ) ;

      iy[i][k] = floor( double(k)/sqrt( double( Nk[i]) ) ) ;
      Y[i][k] = iy[i][k]/sqrt( double (Nk[i]) ) ;
    }

  cout << "ix " ;
  for(int k=0;k<10;k++) 
    cout << ix[0][k] << " " ;
  cout << endl ;
  for(int k=0;k<10;k++) 
    cout << ix[0][Nk[0]-1-k] << " " ;
  cout << endl ;

  cout << "iy " ;
  for(int k=0;k<10;k++) 
    cout << iy[0][k] << " " ;
  cout << endl ;
  for(int k=0;k<10;k++) 
    cout << iy[0][Nk[0]-1-k] << " " ;
  cout << endl ;

  cout << "X " ;
  for(int k=0;k<10;k++) 
    cout << X[0][k] << " " ;
  cout << endl ;
  for(int k=0;k<10;k++) 
    cout << X[0][Nk[0]-1-k] << " " ;
  cout << endl ;

  cout << "Y " ;
  for(int k=0;k<10;k++) 
    cout << Y[0][k] << " " ;
  cout << endl ;
  for(int k=0;k<10;k++) 
    cout << Y[0][Nk[0]-1-k] << " " ;
  cout << endl ;
  
  double p=1. ;
  bool BASELINE=0 ;
  if(abs(Iext[ndI]-m0*sqrt(K)*BL)<=0.0001) {
    cout << "Baseline => No Input Modulation : " << endl ;
    BASELINE=1 ;
  }

  vector<double> IextBL(nbpop) ;
  if(nbpop==1)
    IextBL[0] = BL*sqrt(K)*m0 ;
  else{
    IextBL[0] = Iext[0] ;
    IextBL[1] = BL*sqrt(K)*m0 ;
  }

  p = Iext[ndI] - IextBL[ndI] ;

  cout << "Baseline " ;
  for(int i=0;i<nbpop;i++)
    cout << IextBL[i]/sqrt(K)/m0 << " " ;
  cout << endl ;

  if(BASELINE==0)  {
    cout << "Perturbation " ;
    for(int i=0;i<nbpop;i++)
      cout << Iext[i]/sqrt(K)/m0 << " " ;
  }

  random_device rd ;
  default_random_engine gen( rd() ) ;
  uniform_real_distribution<double> unif(0,1) ;

  for(int i=0;i<nbpop;i++) {
    for(int j=0;j<Nk[i];j++) {
      
      /* if(i==ndI && BASELINE==0) // Modulated input only onto I */
      /* 	Jext[i][j] = p*Wrapped_Gaussian(X[i][j]-X[i][Nk[i]/2],Cff,100,L)*Wrapped_Gaussian(Y[i][j]-Y[i][Nk[i]/2],Cff,100,L) ; */
      /* if(i==ndI && BASELINE==0 && abs(2.*M_PI*(X[i][j]-X[i][Nk[i]/2]))<=Cff) // RectIn */
      /* 	Jext[i][j] = IextBL[i] + p ; */
       
      /* if(i==ndI && BASELINE==0 && (X[i][j]-.5) <= 4*Cff/L && Y[i][j]==Y[i][0]  && unif(gen)<=5 ) //Truncated Gaussian */
      /* 	Jext[i][j] = p*Wrapped_Gaussian(X[i][j]-.5,Cff,100,L) ; */
      if(i==ndI && BASELINE==0 && (X[i][j]-.5)*(X[i][j]-.5)/(4.*Cff/L)/(4.*Cff/L) + (Y[i][j]-.5)*(Y[i][j]-.5)/(4.*Cff/L)/(4.*Cff/L)*10.*10. <= 1. && unif(gen)<=5 ) //Truncated Gaussian
      	Jext[i][j] = p*Wrapped_Gaussian(X[i][j]-.5,Cff,100,L)*Wrapped_Gaussian(Y[i][j]-.5,Cff/10.,100,L) ;
      /* if(i==ndI && BASELINE==0 && sqrt( (X[i][j]-.5)*(X[i][j]-.5) + (Y[i][j]-.5)*(Y[i][j]-.5) ) <= 4.*Cff/L && unif(gen)<=5 ) //Truncated Gaussian */
      /* 	Jext[i][j] = p*Wrapped_Gaussian(X[i][j]-.5,Cff,100,L)*Wrapped_Gaussian(Y[i][j]-.5,Cff,100,L) ; */
      // Jext[i][j] = p ;
      else	
	Jext[i][j] = 0. ;
      
    }
  }
  cout << endl ;


  for(int i=0;i<nbpop;i++) {
    for(int j=0;j<10;j++) 
      cout << Jext[i][j] << " " ;
    cout << endl ;
  }

  delete [] ix ;
  delete [] iy ;
  delete [] X ;
  delete [] Y ;
}

//////////////////////////////////////////////////////////////////

#endif
