#ifndef __NETUTILS__
#define __NETUTILS__

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

using namespace:: std ; 

void printProgress (double percentage) {
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush (stdout);
}

void ProgressBar(double progress) {
  int barWidth = 50;
  std::cout << "[";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos)  std::cout << "\u25a0"; 
    else std::cout << " ";
  }
  std::cout << "] " << int(progress * 100.0) << " %\r";
  std::cout.flush();
  if(progress == 1.) std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////
//Functions
///////////////////////////////////////////////////////////////////////

double phi(int k, int N) {
  return ( (double) k ) *M_PI / (double) N ;
}


double moy (vector<double> v) {
  double sum = accumulate(v.begin(), v.end(), 0.0);
  return sum / v.size();
}

double var (vector<double> v) {
  double sum = accumulate(v.begin(), v.end(), 0.0);
  double mean = sum / v.size();
  double sq_sum = inner_product(v.begin(), v.end(), v.begin(), 0.0);
  return sq_sum / v.size() - mean * mean ;
}

double CV2 (vector<double> v) {
  vector<double> diff ;
  vector<double> sum ;
  
  for(size_t i=0;i<v.size();i++) {
    diff.push_back(fabs(v[i]-v[i+1])) ;
    sum.push_back(v[i]+v[i+1]) ;
  }
  
  return 2.*moy(diff)/moy(sum) ;
}

///////////////////////////////////////////////////////////////////////

void getIextBL(int nbpop, string dir, double* &IextBL) {

  IextBL = new double [nbpop]() ;
  
  string stream_path =  "../Parameters/"+to_string(nbpop)+"pop/"+ dir + "/Param.txt" ;
  cout << stream_path << endl ;
  ifstream stream(stream_path.c_str(), ios::in);
  string line=" "; 
  getline(stream,line);
  line.erase(0,5) ;
  char *end ;
  
  cout << "IextBL : " ;    
  IextBL[0] = strtod(line.c_str(),&end);
  cout << IextBL[0] << " ";
  
  int i=1 ;
  while(i<nbpop) {
    /* cout << i << " " ; */
    IextBL[i] = strtod(end,&end);
    cout << IextBL[i] << " ";
    i++ ;
  }    
  cout << endl;
 
}

///////////////////////////////////////////////////////////////////////

void Set_Parameters(int argc , char** argv, string &dir, int &nbpop, int &N, double &K, double &g, double* &Iext, double* &IextBL) {

  if(argv[1] != NULL) {
    
    dir = argv[1] ;
    nbpop = (int) atoi(argv[2]) ;
    N = (int) atoi(argv[3]) ;
    K = (double) atof(argv[4]) ;
    g = (double) atof(argv[5]) ;

    getIextBL(nbpop,dir,IextBL) ;
    Iext = new double [nbpop] ;      
 
    cout << "Iext : " ;
    if(!argIext) {
      for(int i=0;i<nbpop;i++) {
	Iext[i] = IextBL[i] ;
	cout << Iext[i] << " " ;
      }
      cout << endl ;
    }
    else {
      for(int i=0;i<nbpop;i++) {
	Iext[i] = (double) atof(argv[i+6]) ;
	cout << Iext[i] << " " ;
      }
      cout << endl ;    
    }

  }

  else {
    cout << "Directory ? " ;
    cin >> dir ;
    cout << "nbpop ? " ;
    cin >> nbpop ;
    cout << "total number of neurons ? " ;
    cin >> N ;
    cout << "K ? " ;
    cin >> K ;
    cout << "gain g ? " ;
    cin >> g ;
    cout << "External Inputs ? " ;
    Iext = new double [nbpop] ;      
    for(int i=0;i<nbpop;i++) 
      cin >> Iext[i] ;
  }

}

///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////

void Set_Parameters_Ka(int argc , char** argv, string &dir, int &nbpop, int &N, double* &K, double &g, double* &Iext) {

  if(argv[1] != NULL) {
    dir = argv[1] ;
    nbpop = (int) atoi(argv[2]) ; // total number of neurons prefactor
    N = (int) atoi(argv[3]) ;

    K = new double[nbpop] ;
    for(int i=0;i<nbpop;i++)       
      K[i] = (double) atof(argv[i+4]) ;

    g = (double) atof(argv[5+nbpop]) ;
    
    Iext = new double [nbpop] ;      
    for(int i=0;i<nbpop;i++) 
      Iext[i] = (double) atof(argv[i+nbpop+6]) ;
  }
}

///////////////////////////////////////////////////////////////////////

void Import_Connectivity_Parameters(int nbpop,double** &J,string dir) {

  J = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    J[i] = new double[nbpop] ;
      
  cout << "Reading Connectivity from : " ;
  string Jparam = "../Parameters/" + to_string(nbpop)+"pop/"+ dir +"/Jparam.txt" ;
  cout << Jparam << endl;

  struct stat buffer;   
  if (stat (Jparam.c_str(), &buffer) == 0) {
    FILE *Jfile ;
    Jfile = fopen(Jparam.c_str(),"r");
      
    int dum = 0 ;
    for(int i=0;i<nbpop;i++) 
      for(int j=0;j<nbpop;j++) 
	dum = fscanf(Jfile, "%lf", &J[i][j]) ; 
    
    fclose(Jfile) ;
  }
  else {
    cout << "ERROR ... Jparam.txt not found" << endl ;
    exit(-1) ;
  }

}

///////////////////////////////////////////////////////////////////////

void Import_Synaptic_Parameters(int nbpop,double** &Tsyn,string dir) {

  Tsyn = new double*[nbpop] ;      
  for(int i=0;i<nbpop;i++)
    Tsyn[i] = new double[nbpop] ;
      
  cout << "Reading Synaptic Cst from : " ;
  string Tparam = "../Parameters/" + to_string(nbpop)+"pop/"+ dir +"/Tparam.txt" ;
  cout << Tparam << endl;

  struct stat buffer;   
  if (stat (Tparam.c_str(), &buffer) == 0) {
 
    FILE *Tfile ;
    Tfile = fopen(Tparam.c_str(),"r");
      
    int dum = 0 ;
    for(int i=0;i<nbpop;i++) 
      for(int j=0;j<nbpop;j++) 
	dum = fscanf(Tfile, "%lf", &Tsyn[i][j]) ; 

    fclose(Tfile) ;
  }
  else {
    cout << "Tparam.txt not found" << endl ;
    for(int i=0;i<nbpop;i++) 
      for(int j=0;j<nbpop;j++) {
	if(j==0)
	  Tsyn[i][j] = 3. ;
    	if(j!=0)
	  Tsyn[i][j] = 2. ;    
      }
  }
  
  for(int i=0;i<nbpop;i++) {
    for(int j=0;j<nbpop;j++) 
      cout << Tsyn[i][j] << " " ;
    cout << endl ;
  }
  
}

///////////////////////////////////////////////////////////////////////

void Import_Inputs_number(int nbpop, double K, vector<vector<double> > &Kk, string dir) {
      
  cout << "Reading Inputs number from : " ;
  string Kparam = "../Parameters/" + to_string(nbpop)+"pop/"+ dir +"/Kparam.txt" ;
  cout << Kparam << endl;

  struct stat buffer;   
  if (stat (Kparam.c_str(), &buffer) == 0) {
    FILE *Kfile ;
    Kfile = fopen(Kparam.c_str(),"r");
      
    int dum = 0 ;
    for(int i=0;i<nbpop;i++) 
      for(int j=0;j<nbpop;j++) 
	dum = fscanf(Kfile, "%lf", &Kk[i][j]) ; 
    fclose(Kfile) ;
    for(int i=0;i<nbpop;i++) 
      for(int j=0;j<nbpop;j++)
	Kk[i][j]=Kk[i][j]*K ;
  }
  else {
    cout << "Kparam.txt not found" << endl ;
    for(int i=0;i<nbpop;i++) {
	if(nbpop==2) {
	  Kk[i][0]= 70.*K/100. ;
	  Kk[i][1]= 30.*K/100. ;
	}
	if(nbpop==3) {
	  Kk[i][0]= 70.*K/100. ;
	  Kk[i][1]= 20.*K/100. ;
	  Kk[i][2]= 10.*K/100. ;
	}
	if(nbpop==4) {
	  Kk[i][0]= 70.*K/100. ;
	  Kk[i][1]= 15.*K/100. ;
	  Kk[i][2]= 10.*K/100. ;
	  Kk[i][3]= 5.*K/100. ;
	}
      }
  }
}

///////////////////////////////////////////////////////////////////////

void Import_Connectivity_Matrix(int nbpop,int N,string Jpath,int* &nbPost,unsigned long int* &idxPost,int* &IdPost) {
  
  cout << "Importing Connectivity from : " << endl ;
 
  string nbpath = Jpath ;
  string idxpath = Jpath ;
  string Idpath = Jpath ;
  
  if(IF_Nk) {
    nbpath += "/nbPost_Nk.dat" ;
    idxpath += "/idxPost_Nk.dat" ;
    Idpath += "/IdPost_Nk.dat" ; 
  }
  else {
    nbpath += "/nbPost.dat" ;
    idxpath += "/idxPost.dat" ;
    Idpath += "/IdPost.dat" ; 
  }

  cout << nbpath << endl ;
  cout << idxpath << endl ;
  cout << Idpath << endl ;

  struct stat buffer;   

  if (stat (Idpath.c_str(), &buffer) == 0 || stat (nbpath.c_str(), &buffer) == 0 || stat (idxpath.c_str(), &buffer) == 0) {
    
    N=N*10000 ;
    
    nbPost = new int [N] ;
    idxPost = new unsigned long int [N] ;

    FILE *fnbPost, *fidxPost, *fIdPost ;
    
    int dum ;

    fnbPost = fopen(nbpath.c_str(), "rb") ;
    dum = fread(&nbPost[0], sizeof nbPost[0], N , fnbPost);  
    fclose(fnbPost);

    fidxPost = fopen(idxpath.c_str(), "rb") ;
    dum = fread(&idxPost[0], sizeof idxPost[0], N , fidxPost);
    fclose(fidxPost);
    
    unsigned long int nbposttot = 0 ;
    for(int j=0 ; j<N; j++)
      nbposttot += nbPost[j] ;

    IdPost = new int [nbposttot] ;
    
    fIdPost = fopen(Idpath.c_str(), "rb");
    dum = fread(&IdPost[0], sizeof IdPost[0], nbposttot , fIdPost); 
    fclose(fIdPost);
  }
  else {
    cout << "ERROR : nbPost.dat or idxPost.dat or IdPost.dat not found" << endl ;
    exit(-1) ;
  }

}

///////////////////////////////////////////////////////////////////////

void Check_Connectivity_Matrix(string Jpath,int nbpop, int N, double K) {

  vector<string> List(4) ;
  List[0] = 'E' ;
  List[1] = 'I' ;
  List[2] = 'S' ;
  List[3] = 'V' ;

  string sCon = Jpath ;
  double nbPres = 0.;
  int tot = 0 ;
  double val;

  for(int i=0;i<nbpop;i++) {	
    for(int j=0;j<nbpop;j++) {
      if(IF_Nk)
	sCon = Jpath + "/Con"+List[i]+List[j]+"_Nk.txt" ; 
      else
	sCon = Jpath + "/Con"+List[i]+List[j]+".txt" ; 

      ifstream fCon(sCon.c_str(), ios::in);

      nbPres = 0. ;
      tot = 0 ;
      while(fCon>>val ) {
	nbPres +=val ;
	tot += 1 ;
      }

      string nbPresij = "nbPres"+List[i]+List[j] ;
      cout << nbPresij << "=" << nbPres/(double)tot << "\t" ;

      sCon.clear() ;
      fCon.close();
    }
    cout << endl;
  }

}

///////////////////////////////////////////////////////////////////////

void Save_Parameters(string dir, int nbpop, int N, int *Nk, double K, double** J, double *Iext, const double *Tm, double **Tsyn, string path) {

  string fparam = path + "/Param.txt" ;
  ofstream param(fparam.c_str(), ios::out | ios::ate);

  param << "Number of neurons " ;
  param << N ;
  for(int i=0;i<nbpop;i++)
    param << Nk[i] << " " ;
  param << endl;

  param << "Number of Inputs K " ;
  param << K << endl ;

  param << "dt=" << dt << " "<< "duration=" << duration << " " ;
  param << "Tst=" << Tst << " " << "Tw=" << Tw << endl ; 

  param << "External Inputs" << endl ;
  param << "Iext " << endl ;
  for(int i=0;i<nbpop;i++)
    param << Iext[i] << " " ;
  param << endl ;
 
  param << "Synaptic Strength" << endl ;
  for(int i=0;i<nbpop;i++) {
    for(int j=0;j<nbpop;j++) 
      param << J[i][j] << " " ;
    param << endl;
  }  
      
  param << "Membrane Time Constantes" << endl ;
  for(int i=0;i<nbpop;i++) 
    param << Tm[i] << " " ;
  param << endl ; 

  param << "Synaptic Time Constantes" << endl ;
  for(int i=0;i<nbpop;i++) {
    for(int j=0;j<nbpop;j++) 
      param << Tsyn[i][j]  << " " ;
    param << endl;
  }  
  
  param.close();
}

///////////////////////////////////////////////////////////////////////

void Save_Parameters_Ka(string dir, int nbpop, int N, int *Nk, double *K, double** J, double *Iext, const double *Tm, double **Tsyn, string path) {

  string fparam = path + "/Param.txt" ;
  ofstream param(fparam.c_str(), ios::out | ios::ate);

  vector<char> List(4) ;
  List[0] = 'E' ;
  List[1] = 'I' ;
  List[2] = 'S' ;
  List[3] = 'V' ;

  param << "Number of neurons " ;
  param << N ;
  for(int i=0;i<nbpop;i++)
    param <<" "<< "N"<< List[i] <<"="<< Nk[i] ;
  param << endl;

  param << "Number of Inputs K " ;
  for(int i=0;i<nbpop;i++)
    param << "K" << List[i] << K[i] << " "; 
  param << endl ;

  param << "dt=" << dt << " "<< "duration=" << duration << " " ;
  param << "Tst=" << Tst << " " << "Tw=" << Tw << endl ; 

  param << "External Inputs" << endl ;
  for(int i=0;i<nbpop;i++)
    param << List[i]<<"="<< Iext[i] << " " ;
  param << endl ;
 
  param << "Synaptic Strength" << endl ;
  for(int i=0;i<nbpop;i++) {
    param << "//" ;
    for(int j=0;j<nbpop;j++) 
      param << J[i][j] << " " ;
    param << "//" << endl;
  }  
      
  param << "Membrane Time Constantes" << endl ;
  for(int i=0;i<nbpop;i++) 
    param << Tm[i] << " " ;
  param << endl ; 

  param << "Synaptic Time Constantes" << endl ;
  for(int i=0;i<nbpop;i++) {
    param << "//" ;    
    for(int j=0;j<nbpop;j++) 
      param << Tsyn[i][j]  << " " ;
    param << "//" << endl;
  }  
  
  param.close();
}

///////////////////////////////////////////////////////////////////////

void CreateDir(string dir, int nbpop, int N, double K, double g, string &path) {

  string mkdirp = "mkdir -p " ;
  path += "Simulations/"+ to_string(nbpop)+"pop/" + dir + "/N" + to_string(N) ;

  char cK[10] ;
  sprintf(cK,"%0.0f",K) ;
  string sK = string(cK) ;

  char cg[10] ;
  sprintf(cg,"%0.2f",g) ;
  string sg = string(cg) ;

  path += "/K"+ sK + "/g" + sg ;  
  mkdirp += path ;

  const char * cmd = mkdirp.c_str();
  const int dir_err = system(cmd);

  if(-1 == dir_err) 
    cout << "error creating directories" << endl ;

  cout << "Created directory : " ;
  cout << path << endl ;

}

void CreateDir_Ka(string dir, int nbpop, int N, double *K, double g, string &path) {

  string mkdirp = "mkdir -p " ;
  path += "Simulations/"+to_string(nbpop)+"pop/"+ dir + "/N" + to_string(N) + "/" ;

  string popList[4] = {"E","I","S","V"} ;  
  if(nbpop==1) 
    popList[0] = "I" ;

  char cK[10] ;
  string sK ;
  for(int i=0;i<nbpop;i++) {
    sprintf(cK,"%0.0f",K[i]) ;
    sK = string(cK) ;
    path += "K"+ popList[i] + sK ;
  }

  char cg[10] ;
  sprintf(cg,"%0.2f",g) ;
  string sg = string(cg) ;

  path += "/g" + sg ;  
  mkdirp += path ;

  const char * cmd = mkdirp.c_str();
  const int dir_err = system(cmd);

  if(-1 == dir_err) 
    cout << "error creating directories" << endl ;

  cout << "Created directory : " ;
  cout << path << endl ;

}

void CreateDir_Iext(int nbpop,int nPrtr, double I, string &path) {

  char cI[10] ;
  sprintf(cI,"%0.4f",I) ;

  string sI = string(cI) ;

  string popList[4] = {"E","I","S","V"} ;
  if(nbpop==1)
    popList[0] = "I" ;

  string spop = popList[nPrtr] ;

  string mkdirp = "mkdir -p " ;
  path = path + "/Prtr_" + spop + "/Iext_" + spop + sI ;
  mkdirp += path ;

  const char * cmd = mkdirp.c_str();
  const int dir_err = system(cmd);

  if(-1 == dir_err) 
    cout << "error creating directories" << endl ;

  cout << "Created directories : " ; 
  cout << path << endl ;

}

///////////////////////////////////////////////////////////////////////

void nbNeurons(int nbpop, int &N, int* &Nk) {
  N = N*10000 ;
  /* N = N*2500 ; */

  Nk = new int [nbpop]() ;
  for(int i=0;i<nbpop;i++)
    Nk[i] = N/nbpop ;

  if(IF_Nk) {
    if(nbpop==2) {
      Nk[0]= int(N*80./100.) ;
      Nk[1]= int(N*20./100.) ;
    }
    if(nbpop==3) {
      Nk[0]= int(N*70./100.) ;
      Nk[1]= int(N*15/100.) ;
      Nk[2]= int(N*15./100.) ;

      /* Nk[0]= int(N/3.) ; */
      /* Nk[1]= int(N/3.*25./100.) ; */
      /* Nk[2]= int(N/3.*75./100.) ; */

      /* Nk[0]= int(N/3.) ; */
      /* Nk[1]= int(N/3.*90./100.) ; */
      /* Nk[2]= int(N/3.*10./100.) ; */
    }
    if(nbpop==4) {
      Nk[0]= int(N*70./100.) ; 
      Nk[1]= int(N*10./100.) ; 
      Nk[2]= int(N*10./100.) ; 
      Nk[3]= int(N*10./100.) ; 
    }
  }

}

void cptNeurons(int nbpop, int *Nk, int *&Cpt) {
  Cpt = new int [nbpop+1]() ;
  
  for(int i=0;i<nbpop+1;i++) { 
    for(int j=0;j<i;j++) 
      Cpt[i] += Nk[j] ; 
    cout <<"Cpt="<< Cpt[i] << " ";
  }
  cout << endl ;
}

///////////////////////////////////////////////////////////////////////

void Import_Ring_Orientations(int nbpop, int N, double K, vector<double> &phi){

  string Jpath = "../" ;

  char cK[10] ;
  sprintf(cK,"%0.0f",K) ;
  string sK = string(cK) ;

  Jpath += to_string(nbpop)+"pop/Connectivity/N"+to_string(N/10000)+"/K"+ sK ;

  cout << "Importing Orientations from :" ;

  string str_phi = Jpath + "/ring.dat" ;
  cout << str_phi << endl ;

  struct stat buffer;   
  if (stat (str_phi.c_str(), &buffer) == 0) {
    FILE *fphi ;
    fphi = fopen(str_phi.c_str(), "rb");
    
    int dum = 0 ;
    dum = fread(&phi[0], sizeof phi[0], phi.size(), fphi);  
    
    fclose(fphi);
  }
  else
    cout << "Ring.dat not found ..." << endl ;
}

#endif
