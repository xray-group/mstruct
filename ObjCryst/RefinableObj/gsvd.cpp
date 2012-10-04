// GPL - license
// C-LAPACK - BSD license
// (c) Zdenek


// g++ -I/usr/include/newmat gsvd.cpp -lnewmat -llapack -lblas

// ----------------------- GSVD defines ---

//#define GSVD_TEST_MAIN
//#define GSVD_INTERNAL_INFO

// ----------------------- Newmat ---
#define WANT_STREAM                  // include.h will get stream fns
#define WANT_MATH                    // include.h will get math fns
                                     // newmatap.h will get include.h

#include "newmat/newmatap.h"                // need matrix applications

#include "newmat/newmatio.h"                // need matrix output routines

#ifdef use_namespace
using namespace NEWMAT;              // access NEWMAT namespace
#endif

// ----------------------- LAPACK DGGSVD ---
//#include "clapack.h"               // Some subroutines declaration are missing in clapack.h
                                     // file distributed with ATLAS in Debian 
#define integer         int
#define doublereal      double

// man DGGSVD
extern "C" int dggsvd_(char *jobu, char *jobv, char *jobq, integer *m, 
		       integer *n, integer *p, integer *k, integer *l, doublereal *a, 
		       integer *lda, doublereal *b, integer *ldb, doublereal *alpha, 
		       doublereal *beta, doublereal *u, integer *ldu, doublereal *v, integer 
		       *ldv, doublereal *q, integer *ldq, doublereal *work, integer *iwork, 
		       integer *info);

#ifdef _MSC_VER

#ifdef __USE_FAKED_GSVD__

/* dummy dggsvd_(...) function */
int dggsvd_(char *jobu, char *jobv, char *jobq, integer *m, 
		       integer *n, integer *p, integer *k, integer *l, doublereal *a, 
		       integer *lda, doublereal *b, integer *ldb, doublereal *alpha, 
		       doublereal *beta, doublereal *u, integer *ldu, doublereal *v, integer 
		       *ldv, doublereal *q, integer *ldq, doublereal *work, integer *iwork, 
		       integer *info)
{
	return -1;
}

#endif /* __USE_FAKED_GSVD__ */

#endif /* _MSC_VER */

// ----------------------- other general C++ ---
// using std::vector for automatic memory allocations for arrays   
#include <vector>

using std::vector;

// ----------------------- GSVD ---

/* Generalised SVD decomposition of matrices A and B.
 *
 *   A = U * Sigma1 * [0, R] * Q'
 *   B = V * Sigma2 * [0, R] * Q'
 *
 *   eR = [0, R] ;
 *   r = rank([A; B]) ; k = rank(B) ; l = r-k ;
 */
void GSVD(const Matrix &A, const Matrix &B,
	  Matrix &Sigma1, Matrix &Sigma2,
	  Matrix &U, Matrix &V, Matrix &eR, Matrix &Q,
	  int &K, int &L, vector<int> &ind)
{
  int M = A.Nrows();
  int N = A.Ncols();
  int P = B.Nrows();
  // TODO:: Check B.Ncols()==N
  K = 0, L = 0;
  int LDA = M, LDB = P;

  #ifdef GSVD_INTERNAL_INFO
  {
    cout<<"\n"<<"  --- GSVD decomposition of matrices A and B ---\n";
    cout<<"\t--- A ("<<A.Nrows()<<"x"<<A.Ncols()<<") ---\n";
    cout<<setprecision(4)<<scientific<<showpoint;
    cout<<setw(12) << A;
    cout<<"\t\tMaxAbs: " << A.MaximumAbsoluteValue() << "\n";
    cout<<"\t--- B ("<<B.Nrows()<<"x"<<B.Ncols()<<") ---\n";
    cout<<setprecision(4)<<scientific<<showpoint;
    cout<<setw(12) << B;
    cout<<"\t\tMaxAbs: " << B.MaximumAbsoluteValue() << "\n";
  }
  #endif // GSVD_INTERNAL_INFO

  // Space for output matrices
  vector<double> vAlpha(N,0.); // in Fortran should be 1D array
  vector<double> vBeta(N,0.);

  // We must transpose data in A and B, becouse of column-major
  // data order in Fortran GSVD
  Matrix AR = A.t();
  Matrix BR = B.t();
  if(P<1) { BR = Matrix(1,N); }

  char JOBU = 'U', JOBV = 'V', JOBQ = 'Q';
  
  int LDU = max(1,M), LDV = max(1,P), LDQ = max(1,N);
  U.ReSize(LDU,M);
  V.ReSize(LDV,P);
  Q.ReSize(LDQ,N);

  int wk = max(3*N,M); wk = max(wk,P); wk += N;
  vector<double> WORK(wk,0.);

  vector<int> & IWORK = ind; // alias
  IWORK = vector<int>(N,0);

  int INFO = 0;
  int RESULT = 0;

  RESULT = dggsvd_(&JOBU, &JOBV, &JOBQ, &M, &N, &P, &K, &L,
		   AR.Store(), &LDA, BR.Store(), &LDB,
		   &(vAlpha[0]), &(vBeta[0]), U.Store(), &LDU, V.Store(), &LDV,
		   Q.Store(), &LDQ, &(WORK[0]), &(IWORK[0]), &INFO);

  #ifdef GSVD_INTERNAL_INFO
  {
    cout<<"\t--- GSVD calculated ---\n";
    cout<<"\tk: "<<K<<", l: "<<L<<", info: "<<INFO<<", result: "<<RESULT<<"\n";
  }
  #endif // GSVD_INTERNAL_INFO

  // assign matrices Sigma1 and Sigma2
  Sigma1.ReSize(M,K+L); Sigma1 = 0.;
  Sigma2.ReSize(P,K+L); Sigma2 = 0.;
  for(int i=1; i<=K; i++) Sigma1(i,i) = 1.;
  if(M-K-L>=0) {       
    for(int i=(K+1); i<=(K+L); i++) {
      Sigma1(i,i) = vAlpha[i-1];
      Sigma2(i-K,i) = vBeta[i-1];
    }       
  } else {
    for(int i=(K+1); i<=M; i++) {
      Sigma1(i,i) = vAlpha[i-1];
      Sigma2(i-K,i) = vBeta[i-1];
    }
    for(int i=(M+1); i<=(K+L); i++)
      Sigma2(i-K,i) = vBeta[i-1];
  }

  #ifdef GSVD_INTERNAL_INFO
  {
    cout<<"\t--- Sigma1 ("<<Sigma1.Nrows()<<"x"<<Sigma1.Ncols()<<") ---\n";
    cout<<setprecision(4)<<scientific<<showpoint;
    cout<<setw(12) << Sigma1;
    cout<<"\t--- Sigma2 ("<<Sigma2.Nrows()<<"x"<<Sigma2.Ncols()<<") ---\n";
    cout<<setprecision(4)<<scientific<<showpoint;
    cout<<setw(12) << Sigma2;
  }
  #endif // GSVD_INTERNAL_INFO
  
  // transpose U, V, Q matrices from the column-major Fortran format
  U = U.t();
  V = V.t();
  Q = Q.t();
  // transpose also AR and BR matrices from the column-major Fortran format
  AR = AR.t();
  BR = BR.t();

  // extract and assambly eR (extended R) matrix
  eR.ReSize(K+L,N); eR = 0.;
  if(M-K-L>=0) {
    // copy data from AR matrix
    for(int i=1; i<=(K+L); i++)
      //for(int j=N-K-L+1; j<=N)
      for(int j=(N-K-L+i); j<=N; j++) eR(i,j) = AR(i,j);
  } else {
    // copy data from AR and BR matrices
    for(int i=1; i<=M; i++)
      //for(int j=N-K-L+1; j<=N)
      for(int j=(N-K-L+i); j<=N; j++) eR(i,j) = AR(i,j);
    for(int i=(M-K+1); i<=L; i++)
      //for(int j=N-K-L+1; j<=N)
      for(int j=(N+M-K-L+1); j<=N; j++) eR(i,j) = BR(i,j);
  }

  #ifdef GSVD_INTERNAL_INFO
  {
    cout<<"\t--- U ("<<U.Nrows()<<"x"<<U.Ncols()<<") ---\n";
    cout<<setprecision(4)<<scientific<<showpoint;
    cout<<setw(12) << U;
    
    cout<<"\t--- V ("<<V.Nrows()<<"x"<<V.Ncols()<<") ---\n";
    cout<<setprecision(4)<<scientific<<showpoint;
    cout<<setw(12) << V;

    cout<<"\t--- eR = [0, R] ("<<eR.Nrows()<<"x"<<eR.Ncols()<<") ---\n";
    cout<<setprecision(4)<<scientific<<showpoint;
    cout<<setw(12) << eR;
    
    cout<<"\t--- Q ("<<Q.Nrows()<<"x"<<Q.Ncols()<<") ---\n";
    cout<<setprecision(4)<<scientific<<showpoint;
    cout<<setw(12) << Q;

    cout<<"\t--- ind ("<<ind.size()<<") ---\n";
    cout<<setprecision(4)<<scientific<<showpoint;
    for(int i=0; i<ind.size(); i++) cout<<setw(8) << ind[i];
    cout<<"\n";
  }
  #endif // GSVD_INTERNAL_INFO
   
  #ifdef GSVD_INTERNAL_INFO
  {
    cout<<"\t--- GSVD test: A - U*Sigma1*eR*Q' ---\n";
    Matrix diff = A - U*Sigma1*eR*Q.t();
    cout<<setprecision(3)<<scientific<<showpoint;
    cout<<setw(10) << diff;
    cout<<setprecision(3)<<scientific<<showpoint;  
    cout<<"\t\tMaxAbs: " << diff.MaximumAbsoluteValue() << "\n";

    cout<<"\t--- GSVD test: B - V*Sigma2*eR*Q' ---\n";
    diff = B - V*Sigma2*eR*Q.t();
    cout<<setprecision(3)<<scientific<<showpoint;
    cout<<setw(10) << diff;
    cout<<setprecision(3)<<scientific<<showpoint;  
    cout<<"\t\tMaxAbs: " << diff.MaximumAbsoluteValue() << "\n";
    
    cout<<flush;
  }
  #endif // GSVD_INTERNAL_INFO
}

// ----------------------- main for tests ---

#ifdef GSVD_TEST_MAIN

void test1 (Matrix &A, Matrix &B,
	    Matrix &Sig1, Matrix &Sig2,
	    Matrix &U, Matrix &V, Matrix &eR, Matrix &Q,
	    int &k, int &l, vector<int> &ind)
{
  cout<<"\t--- TEST 1 ---\n";
  
  A.ReSize(3,2);
  A(1,1) = 1.; A(1,2) = 0.;
  A(2,1) = 1.; A(2,2) = 1.;
  A(3,1) = 1.; A(3,2) = 3.;
  B.ReSize(1,2);
  B(1,1) = 1.; B(1,2) = 4.;
  
  GSVD(A,B,Sig1,Sig2,U,V,eR,Q,k,l,ind);
  
  if( (k+l)!=A.Ncols() ) {
    cerr << "Math error: rank([A; B]) neq nb.ofColms A."<<endl;
    exit(0);
  }

  Matrix W = Q * eR.i();

  ColumnVector b(3);
  b(1) = 0.; b(2) = 8.; b(3) = 8.;
  ColumnVector d(1);
  d(1) = 20.;

  ColumnVector x(2); x = 0.;

  for(int i=1; i<=k; i++) {
    ColumnVector u = U.Column(i);
    x += DotProduct(u.t(),b)/Sig1(i,i) * W.Column(i);
  }

  for(int i=(k+1); i<=(k+l); i++) {
    ColumnVector v = V.Column(i-k);
    x += DotProduct(v.t(),d)/Sig2(i-k,i) * W.Column(i);
  }
  
  cout << "\t--- Constrained solution x ---\n";
  cout << scientific << showpoint << setprecision(4);
  cout << " " << setw(10) << x.t();
  
  ColumnVector diff = d - B * x;
  
  cout << "\t--- Constrained solution x: (d - B*x)' ---\n";
  cout << scientific << showpoint << setprecision(3);
  cout << " " << setw(10) << diff.t();
  cout << "\t\tMaxAbs: " << diff.MaximumAbsoluteValue() << "\n";
}

int main(int argc, char* argv[])
{
  { // Check if C, Newmat and Fortran data types are of the same size
    cout << "double:     " << sizeof(double) << "\n";
    cout << "float:      " << sizeof(float) << "\n";
    cout << "Real:       " << sizeof(Real) << "\n";
    cout << "doublereal: " << sizeof(doublereal) << "\n";
    cout << "int:        " << sizeof(int) << "\n";
    cout << "integer:    " << sizeof(integer) << "\n";

    if( sizeof(double)!=sizeof(Real) ) {
      cerr << "Compatibility error: Newmat Real is not of double type!" << endl;
      exit(0);
    }
  }

  { // Check low level access to the Newmat matrices
    int n = 5;
    Matrix A(n,n);

    for(int i=1; i<=n; i++)
      for(int j=1; j<=n; j++) A(i,j) = (i-1)*n+j;
  
    cout << "      --- A ---" << "\n";
    cout << setprecision(1) << fixed << showpoint;
    cout << setw(8) << A << endl;

    cout << "      --- A --- low level access" << "\n";
    cout << setprecision(1) << fixed << showpoint;
    Real *p = A.Store();
    for(int i=0; i<n; i++) {
      cout << setw(8);
      for(int j=0; j<n; j++){ cout << *p; p++; cout << setw(9); }
      cout << "\n";
    }

    cout << endl;
  }

  { // test1
    Matrix A, B, Sig1, Sig2, U, V, eR, Q;
    int k, l;
    vector<int> ind;

    test1(A,B,Sig1,Sig2,U,V,eR,Q,k,l,ind);
  }

  return 0;
}

#endif // GSVD_TEST_MAIN
 
