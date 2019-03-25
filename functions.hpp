#ifndef _FUNCTIONS_HPP_DEFINED_
#define _FUNCTIONS_HPP_DEFINED_

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <lapacke.h>


extern int size;
extern double t;
extern double G;

int size;
double t=1;
double G;

using std::cout;
using std::endl;
using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using std::vector;

inline int xc(int i, int sigma=1){return (sigma==1)?std::floor(i/size):std::floor((i-size*size)/size);}
inline int yc(int j, int sigma=1){return (sigma==1)?j%size:(j-size*size)%size;}
inline int periodic(int a, int b, int lim) {int rem = (a+b)%lim; if(rem>=0) return rem; else return rem+lim;}

Eigen::MatrixXd construct_h0_2d(void)
{
  MatrixXd Mc = MatrixXd::Zero(2*size*size, 2*size*size);
  for(int i=0; i<size*size; i++)
  {
    for(int j=0; j<size*size; j++)
    {
      int ix = xc(i); int iy = yc(i); int jx = xc(j); int jy = yc(j);
      if((ix == periodic(jx,1,size)|| ix == periodic(jx,-1,size)) && iy == jy ) Mc(i,j) = -t;
			if(ix == jx && (iy == periodic(jy,1,size) || iy == periodic(jy,-1,size))) Mc(i,j) = -t;				
    }
  }

  for(int i=size*size; i<2*size*size; i++)
  {
    for(int j=size*size; j<2*size*size; j++)
    {
      int ix = xc(i,-1); int iy = yc(i,-1); int jx = xc(j,-1); int jy = yc(j,-1);
      if((ix == periodic(jx,1,size)|| ix == periodic(jx,-1,size)) && iy == jy ) Mc(i,j) = -t;
			if(ix == jx && (iy == periodic(jy,1,size) || iy == periodic(jy,-1,size))) Mc(i,j) = -t;				
    }
  }
  return Mc;
}

MatrixXd construct_h0(void)
{
  MatrixXd Mc = MatrixXd::Zero(2*size,2*size);
  for(int row=0; row <2*size-1; row++) Mc(row,row+1)=Mc(row+1,row)=-t;
  Mc(size-1,0)=Mc(0,size-1)=-t; //PBC
  Mc(2*size-1, size)=Mc(size,2*size-1)= -t; //PBC
  Mc(size,size-1)= Mc(size-1,size)=0;
  return Mc;
}

bool zheev_cpp(MatrixXcd& A, vector<double>& lambda, char eigenvec_choice='N')
{
  int N = A.cols();
  int LDA = A.outerStride();
  int INFO = 0;
  double* w = new  double [N];
  char Nchar = 'N';
  char Vchar = 'V';
  char Uchar = 'U';
  int LWORK = int(A.size())*4;
  __complex__ double* WORK= new __complex__ double [LWORK];
  double* RWORK = new double [3*LDA];

  zheev_( &eigenvec_choice, &Uchar, &N, reinterpret_cast <__complex__ double*> (A.data()), &LDA, w, WORK, &LWORK, RWORK, &INFO );

  lambda.clear();
  for(int i=0; i<N; i++) lambda.push_back(w[i]);

  delete[] w; delete[] RWORK; delete[] WORK;
  return INFO==0;
}

bool zgeev_cpp(MatrixXcd& A, vector<double>& lambda, char eigenvec_choice='N')
{  
  int N = A.cols();
  int LDA = A.outerStride();
  int INFO = 0;
  __complex__ double* w = new __complex__ double [N];
  __complex__ double* vl;
  __complex__ double* vr;
  char Nchar = 'N';
  char Vchar = 'V';
  char Uchar = 'U';
  int LWORK = pow(2, N);
  __complex__ double* WORK= new __complex__ double [LWORK];
  double* RWORK = new double [LWORK];
  
  zgeev_(&Nchar, &eigenvec_choice, &N, reinterpret_cast <__complex__ double*> (A.data()), &LDA, w, vl, &LDA, vr, &LDA, WORK, &LWORK, RWORK, &INFO );

  lambda.clear();
  for(int i=0; i<N; i++) lambda.push_back(__real__ w[i]);

  delete[] w; delete[] RWORK; delete[] WORK;
  return INFO==0;
}

vector <double> stdEigenvalues(MatrixXcd A, bool (*diagonalization_routine)(MatrixXcd&, vector <double>&, char)=&zheev_cpp)
{
  std::vector<double> lambda; 
  if(diagonalization_routine(A,lambda,'N')) return lambda;
}

VectorXd Eigenvalues(MatrixXcd A, bool (*diagonalization_routine)(MatrixXcd&, vector <double>&, char)=&zheev_cpp)
{
  std::vector<double> lambda;
  if(diagonalization_routine(A,lambda,'N'))
 	{
		Eigen::Map<Eigen::ArrayXd> b(lambda.data(),lambda.size());
  	return b;
	}
}

MatrixXcd Eigenvectors(MatrixXcd A, bool (*diagonalization_routine)(MatrixXcd&, vector <double>&, char)=&zheev_cpp)
{
	std::vector<double> lambda;
  if(diagonalization_routine(A,lambda,'V')) return A; 
}

std::pair<MatrixXcd, vector<double> > stdEigenspectrum(MatrixXcd A, bool (*diagonalization_routine)(MatrixXcd&, vector <double>&, char)=&zheev_cpp)
{
  std::vector<double> lambda;
  if(diagonalization_routine(A,lambda,'V')) return make_pair(A,lambda);
}

std::pair<MatrixXcd, VectorXd> Eigenspectrum(MatrixXcd A, bool (*diagonalization_routine)(MatrixXcd&, vector <double>&, char)=&zheev_cpp)
{
  std::vector<double> lambda;
  if(diagonalization_routine(A,lambda,'V'))
 	{
    Eigen::Map<Eigen::ArrayXd> b(lambda.data(),lambda.size());
    return make_pair(A,b);
	}
}

double dos(double omega, const Eigen::VectorXd& eivals)
{
  double density = 0.0;
  for(int i=0; i<eivals.size(); i++)
  {
    density += (G/2)/( pow((G/2),2) + pow((omega-eivals(i)),2) );
  }
  return density/M_PI;
}

#endif