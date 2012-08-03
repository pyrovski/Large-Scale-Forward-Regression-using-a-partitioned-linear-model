#include "pinv.h"
#include "svd.h"

using namespace std;

int pinv(FortranMatrix &M, vector<double> &rhs, FortranMatrix &Mi, 
	 vector<double> &beta, double &tol, int id){
  int n = M.get_n_rows();
  vector<double> S;
  FortranMatrix U, Vt;

  int rank;
  // Initialize SVD components, A = U * S * V^T
  //! @todo this is altering M?
  svd_create(M, U, S, Vt);

  // beta = V' * S^-1 * U' * rhs
  rank = svd_apply(U, S, Vt, /*result=*/beta, rhs);

  if(!id){
    U.writeD("U.dat");
    writeD("S.dat", S);
    Vt.writeD("Vt.dat");
    writeD("beta.dat", beta);
  }

  // Mi = V * S^-1 * Ut
  // S^-1 = 1./S, where S > tol

  /* compute V * 1./S
     S is stored as a vector, but represents an nxn diagonal matrix
  */

  // first compute V from Vt
    
  Vt.transpose_self();

  double maxS = 0.0;
  for(unsigned i = 0; i < n; i++)
    maxS = max(maxS, S[i]); // compute norm(M, 2) = max(S)
  tol = n * numeric_limits<double>::epsilon() * maxS;
  for(unsigned i = 0; i < n; i++)
    if(S[i])
      if(S[i] > tol)
	S[i] = 1.0/S[i];
      else
	S[i] = 0.0;
    
  // emulate matrix-matrix multiply, with second matrix diagonal
  for(unsigned col = 0; col < n; col++)
    cblas_dscal(n,
		S[col],
		&Vt.values[col * n],
		1);

  // compute Mi = V * S^-1 * Ut
  Mi.resize(n, n);
  cblas_dgemm(CblasColMajor,
	      CblasNoTrans, // V_Si is not transposed
	      CblasTrans, // U is transposed
	      n,
	      n,
	      n,
	      1.0,
	      &Vt.values[0],
	      n,
	      &U.values[0],
	      n,
	      0.0,
	      &Mi.values[0],
	      n);

  return rank;
}
