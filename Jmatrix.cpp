#include "librairies.h"
#include "Net_Utils.h"
#include "Matrix_Utils.h"

#define IF_Nk true
#define IF_MATRIX false

using namespace::std ;

clock_t t1=clock();

int main(int argc , char** argv) {
 
  ///////////////////////////////////////////////////////////////////    
  // Parameters 
  ///////////////////////////////////////////////////////////////////    
 
  int nbpop = (int) atoi(argv[1]) ; // number of neuronal populations
  int N = (int) atoi(argv[2]) ; // total number of neurons prefactor, Ntot = N*10000
  double K = (double) atof(argv[3]) ;  // average number of connections

  ///////////////////////////////////////////////////////////////////    
  // Path
  ///////////////////////////////////////////////////////////////////    

  string path = "../" ;
  Create_Path(nbpop,path,N,K) ;

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

  cout <<"Generating Connection Probabilities ..." << endl;
  
  vector<vector<double> >c(nbpop,vector<double> (nbpop)) ;
  for(int i=0;i<nbpop;i++)
    for(int j=0;j<nbpop;j++)
      c[i][j] = K/(double)Nk[j] ; 
  Save_Param(nbpop,path,K,Nk,c) ;
  
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
  
  ///////////////////////////////////////////////////////////////////    

  random_device rd ;
  default_random_engine gen( rd() ) ;

  uniform_real_distribution<double> unif(0,1) ;

  cout << "Check random seed " ;
  for(int i=0;i<10;i++)
    cout << unif(gen) << " " ;
  cout << endl ;

  ///////////////////////////////////////////////////////////////////    
  
  cout <<"Generating vectors nbPost & IdPost ..." <<endl;

  for(int i=0;i<nbpop;i++) 
    for(int k=Cpt[i];k<Cpt[i+1];k++) //Presynaptic neurons
      for(int j=0;j<nbpop;j++)
	for(int l=Cpt[j];l<Cpt[j+1];l++) //Postsynaptic neurons
	  if(unif(gen)<=c[j][i]) { // from uniform to Bernoulli
	    IdPost.push_back(l); 
	    nbPost[k]++ ;
	    nbPreSab[j][i][l-Cpt[j]]++ ;
	  }

  c.clear() ;

  ///////////////////////////////////////////////////////////////////    
  // Average number of Presynaptic neurons
  ///////////////////////////////////////////////////////////////////    
  
  // NumberPres(nbpop,path,Nk,nbPreSab,IF_Nk) ;
  CheckPres(nbpop,path,Nk,nbPreSab) ;
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
