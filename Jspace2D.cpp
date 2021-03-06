#include "Matrix_Utils.h"
#include "Space_Utils.h"
#include "Net_Utils.h"

using namespace::std ;

clock_t t1=clock();

#define L 2.0*M_PI
#define IF_Nk false
#define IF_MATRIX true
#define klim 10 // max bound of the wrapped gaussian sum

int main(int argc , char** argv) {
 
  ///////////////////////////////////////////////////////////////////    
  // Parameters 
  ///////////////////////////////////////////////////////////////////    
 
  int nbpop = (int) atoi(argv[1]) ; // number of neuronal populations
  int N = (int) atoi(argv[2]) ; // total number of neurons prefactor, Ntot = N*10000.
  double K = (double) atof(argv[3]) ;  // average number of connections  
  double *Crec ; // variance of the wrapped
  Crec = new double [nbpop] ;
  for(int i=0;i<nbpop;i++) 
    Crec[i] = (double) atof(argv[4+i]) ;

  ///////////////////////////////////////////////////////////////////    
  // Path
  ///////////////////////////////////////////////////////////////////    

  string path = "../" ;
  Create_Path(nbpop,path,N,K) ;
  CreateDir_SpaceCrec2D(nbpop,path,N,Crec) ;

  ///////////////////////////////////////////////////////////////////    
  // Number of Neurons
  ///////////////////////////////////////////////////////////////////    

  int* Nk ;
  nbNeurons(nbpop,N,Nk,IF_Nk) ;
  
  vector<int> Cpt(nbpop+1) ;
  for(int i=0;i<nbpop+1;i++) {
    for(int j=0;j<i;j++) 
      Cpt[i] += Nk[j] ; 
    cout <<"Cpt="<< Cpt[i] << " ";
  }
  cout << endl ;

  ///////////////////////////////////////////////////////////////////
  // Connection Probabilities
  ///////////////////////////////////////////////////////////////////

  
  // vector<vector<double> >c(nbpop,vector<double> (nbpop)) ;
  // for(int i=0;i<nbpop;i++)
  //   for(int j=0;j<nbpop;j++)
  //     c[i][j] = K/(double)Nk[j] ; 
  // Save_Param(nbpop,path,K,Nk,c) ;

  double **c ;
  Spatial_Connection_Probability_2D(nbpop,N,Nk,Cpt,K,klim,Crec,c,L) ;

  ///////////////////////////////////////////////////////////////////
  // Sparse Vectors 
  ///////////////////////////////////////////////////////////////////

  vector<int> IdPost ; // Id of the post neurons
  vector<int> nbPost(N) ; // number of post neurons
  vector<unsigned long int> idxPost(N) ; // idx of the post neurons
  vector<vector<vector<int> > > nbPreSab(nbpop,vector<vector<int> >(nbpop) ) ;
  
  idxPost[0] = 0 ;

  for(int i=0;i<nbpop;i++)
    for(int j=0;j<nbpop;j++)
      for(int k=0;k<Nk[i];k++)
  	nbPreSab[i][j].push_back(0) ;

  random_device rd ;
  default_random_engine gen( rd() ) ;
  uniform_real_distribution<double> unif(0.,1.) ;

  ///////////////////////////////////////////////////////////////////    

  cout <<"Generating vectors nbPost & IdPost ..." << endl;

  for(int i=0;i<nbpop;i++) 
    for(int k=Cpt[i];k<Cpt[i+1];k++) //Presynaptic neurons i to j
      for(int j=0;j<nbpop;j++)
	for(int l=Cpt[j];l<Cpt[j+1];l++) //Postsynaptic neurons
	  if(unif(gen)<c[l][k]) { // from uniform to Bernoulli, k to l
	    IdPost.push_back(l) ; 
	    nbPost[k]++ ;
	    nbPreSab[j][i][l-Cpt[j]]++ ;
	  }

  delete [] c ;
  ///////////////////////////////////////////////////////////////////    
  // Average number of Presynaptic neurons
  ///////////////////////////////////////////////////////////////////    
  
  NumberPres(nbpop,path,Nk,nbPreSab,IF_Nk) ;
  CheckPres(nbpop,Nk,nbPreSab) ;
  nbPreSab.clear() ;

  ///////////////////////////////////////////////////////////////////    
  // Writing to File
  ///////////////////////////////////////////////////////////////////

  WritetoFile(path,N,IdPost,nbPost,idxPost,IF_Nk) ;

  ///////////////////////////////////////////////////////////////////    
  // Writing Complete Matrix
  ///////////////////////////////////////////////////////////////////

  if(IF_MATRIX)
    WriteMatrix(path,N,IdPost,nbPost,idxPost) ;

  ///////////////////////////////////////////////////////////////////

  IdPost.clear() ;
  idxPost.clear();
  nbPost.clear();

  cout << "Done" << endl ;
}
