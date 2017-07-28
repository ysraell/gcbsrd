
/******************************************************
 *
 *
 *
 * mex -I/usr/include/eigen3 ssimcpp.cpp
 ******************************************************/

#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include "mex.h"


using namespace Eigen;
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  const mwSize *dim;
  dim = mxGetDimensions(prhs[0]);
  
  double *Aptr, *Bptr;
  Aptr = mxGetPr(prhs[0]);
  Bptr = mxGetPr(prhs[1]);
   
  MatrixXd A, B, A2, B2;
  A.resize(dim[0],dim[1]);
  A << Map<MatrixXd>(Aptr, dim[0], dim[1]);
  B.resize(dim[0],dim[1]);
  B << Map<MatrixXd>(Bptr, dim[0], dim[1]);
  
  A2.resize(dim[0],dim[1]);
  B2.resize(dim[0],dim[1]);
    
  const double mA = A.mean();
  const double mB = B.mean();
  
  A2 = (A- mA*MatrixXd::Ones(dim[0], dim[1]));
  B2 = (B- mB*MatrixXd::Ones(dim[0], dim[1]));
  
  const double N = dim[0]*dim[1];
  const double sA = A2.cwiseAbs2().sum()/N;
  const double sB = B2.cwiseAbs2().sum()/N;
  
  double sAB = 0.0;
  uint i, j;
  for (i=0; i<N; i++){
      for (j=0; j<N; j++){
          if (j>i) sAB = sAB+(A(i)-A(j))*(B(i)-B(j));
      }
  }
  sAB = sAB/(N*N);
  
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  double *Out;
  Out = mxGetPr(plhs[0]);
  MatrixXd::Map(Out, 1, 1) << (2*mA*mB+0.001)*(2*sAB+0.009)/((mB*mB+mA*mA+0.001)*(sB+sA+0.009));
  return;
  
}

//EOF