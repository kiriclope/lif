#ifndef __MFRATES__
#define __MFRATES__

#include <iostream>
#include <eigen3/Eigen/Dense>

using namespace Eigen;

void MF_Rates(int nbpop,double* Iext, double** J, double* &RatesMF) {
  
  RatesMF = new double [nbpop]() ;

  MatrixXd A(nbpop,nbpop) ;
  VectorXd b(nbpop) ;
  
  for(int i=0;i<nbpop;i++) {
    b(i) = Iext[i]*m0 ; 
    for(int j=0;j<nbpop;j++) 
      A(i,j) = J[i][j] ; 
  }
  
  /* cout << "Here is the matrix A:\n" << A << endl; */
  /* cout << "Here is the vector b:\n" << b << endl; */
  VectorXd x(nbpop) ; 
  x = A.colPivHouseholderQr().solve(-b) * 1000. ;
  /* cout << "The solution is:\n" << x << endl; */
  
  for(int i=0;i<nbpop;i++)
    RatesMF[i] = x(i) ; 
  
}

#endif
