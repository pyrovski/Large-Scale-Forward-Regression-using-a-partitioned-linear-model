// System header files
#include <sys/time.h>
#include <string>
#include <fstream>
#include <vector>
#include <limits>
#include <algorithm>
#include <cuda.h>

// Local project includes
#include "fortran_matrix.h"
#include "glm.h"
#include "GetPot"
#include "svd.h"
#include "cblas.h"

#include "plm.cu"

using namespace std;


int main()
{
  // Timing variables
  timeval tstart, tstop;

  // Create input file object.  Put the path to your data files here!
  GetPot input_file("reference_glm.in");
  
  // The path on my system to the location of the data files.  Don't forget the trailing
  // slash here, as this will be prepended to the filename below
  string path = "./"; // default path is the current directory.
  path = input_file("path", path.c_str());
  
  // File containing the "population structure".  It is a 4892-by-26 matrix
  string fixed_filename = "fixed.effects.nam.sorted.filtered.dat";

  // File containing the genotypes.  It is a 4892-by-79 matrix.
  string geno_filename = "imputed.marker.chr10.sorted.filtered.dat";

  // File containing the phenotypes.  It is a 4892-by-1 matrix.  The file is
  // arranged in a single column.
  string y_filename = "residuals.chr10.sorted.dat";

  // In Matlab, these sizes are inferred from the data.  In C++, we hard-code them
  // to make reading the data simpler...
  unsigned pop_ind = 4892, fixed_count = 26; // rows, columns of the fixed array
  unsigned geno_ind = 4892, geno_count = 79; // rows, columns of the geno array
  unsigned y_ind = 4892, y_count = 1;        // rows, columns of the y vector

  // Matrix objects for storing the input data
  FortranMatrix fixed(pop_ind, fixed_count);
  FortranMatrix geno(geno_ind, geno_count);
  vector<double> y(y_ind);

  // Begin timing the file IO for all 3 files
  gettimeofday(&tstart, NULL);

  
  // Read the "fixed" array from file
  {
    // Open "fixed" file for reading
    string filename = path + "/" + fixed_filename;
    ifstream fixed_file(filename.c_str());

    if (fixed_file){
      // Loop over all rows and columns, set entries in the fixed matrix
      // double val=99;
      for (unsigned i=0; i<fixed.get_n_rows(); ++i)
	for (unsigned j=0; j<fixed.get_n_cols(); ++j){
	  fixed_file >> fixed(i,j);
	}
    } else {
      cout << "Failed to open file: " << fixed_file << "!!" << endl;
      return 1;
    }
  }

  // Read the geno array from file

  {
    // Open "fixed" file for reading
    ifstream geno_file((path + "/" + geno_filename).c_str());
    
    if (geno_file){
      // Loop over all rows and columns, set entries in the matrix
      for (unsigned i=0; i<geno.get_n_rows(); ++i)
	for (unsigned j=0; j<geno.get_n_cols(); ++j){
	  geno_file >> geno(i,j);
	}
    } else {
      cout << "Failed to open file!!" << endl;
      return 1;
    }
  }
  
  
  // Read the y-array from file.  Currently stored as a vector since that
  // is how it is passed to the glm function, but could be changed to a
  // FortranMatrix with one column...

  {
    // Open "fixed" file for reading
    ifstream y_file((path + "/" + y_filename).c_str());

    if (y_file){
	// Loop over all rows and columns, set entries in the matrix
	for (unsigned i=0; i<y.size(); ++i) {
	  y_file >> y[i];
	}
    } else {
      cout << "Failed to open file!!" << endl;
      return 1;
    }
  }
  
  gettimeofday(&tstop, NULL);
  
  {
    // Compute time taken for IO
    const double io_elapsed_time = 
      (static_cast<double>(tstop.tv_sec  - tstart.tv_sec) +
       static_cast<double>(tstop.tv_usec - tstart.tv_usec)*1.e-6);
    
    cout << "Time required for I/O: " << io_elapsed_time << " s." << endl;
  }
  
  // Version A, Kt a general matrix.
  // Create the Kt matrix.  It has 1 row and fixed_count+2 columns.
  // The entries of Kt are all zero except for the last entry, which is 1.
  //FortranMatrix Kt(1,fixed_count+2);
  //Kt(0, fixed_count+1) = 1.; // Set last entry = 1

  // Version B, Kt assumed a vector.
  vector<double> Kt(fixed_count+2);
   Kt.back() = 1.; // Set last entry = 1
  
  // An array to hold the results of the GLM calculations
  vector<double> Pval(geno_count);
  vector<double> Fval(geno_count);
  vector<double> V2s(geno_count);
  
  // Initialize the X-matrix.  The first column is all ones, the next
  // fixed_count columns are equal to the fixed matrix, and the last
  // column (which changes) is the i'th column of the geno array.
  unsigned n = fixed_count + 1;
  FortranMatrix X(geno_ind/*4892*/, n/*27*/);

  // Fill first column of X with 1's
  for (unsigned i=0; i<X.get_n_rows(); ++i)
    X(i,0) = 1.;

  // Fill next fixed_count columns with the fixed array contents
  memcpy(&X.values[X.get_n_rows()], &fixed.values[0], 
	 fixed_count * X.get_n_rows() * sizeof(double));

  
  // Begin timing the computations
  gettimeofday(&tstart, NULL);

  /*! precompute
    XtX
    XtXi
    rank(XtXi)
    
    for each SNP:
    SNPty
    SNPtSNP
   */
  
  FortranMatrix XtX = matmat(X, X, true, false); 
  XtX.writeD("XtX.dat");

  // Solve (X^T * X)*beta = X^T*y for beta.  Note that X and X^T * X
  // have the same rank.

  // Initialize SVD components, A = U * S * V^T
  vector<double> S;
  FortranMatrix U, Vt, XtXi;

  GLMData glm_data, glm_data_new;
  vector<double> beta;
  vector<double> Xty = matvec(X, y, /*transX=*/true);
  
  writeD("Xty.dat", Xty);

  // Create the SVD of X^T * X 
  svd_create(XtX, U, S, Vt);
  U.writeD("U.dat");
  writeD("S.dat", S);

  Vt.writeD("Vt.dat");
  unsigned rX = svd_apply(U, S, Vt, /*result=*/beta, Xty);

  writeD("beta.dat", beta);

  // XtXi = V * S^-1 * Ut
  // S^-1 = 1./S, where S > tol

  /* compute V * 1./S
     S is stored as a vector, but represents an nxn diagonal matrix
  */

  // first compute V from Vt
    
  Vt.transpose_self();

  double maxS = 0.0;
  for(unsigned i = 0; i < n; i++)
    maxS = max(maxS, S[i]); // compute norm(XtX, 2) = max(S)
  double tol = n * numeric_limits<double>::epsilon() * maxS;
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

  // compute XtXi = V * S^-1 * Ut
  XtXi.resize(n, n);
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
	      &XtXi.values[0],
	      n);

  XtXi.writeD("XtXi.dat");

  // Compute the matrix-vector product, XTy := X' * y.  
  double yty = cblas_ddot(y.size(), &y[0], 1, &y[0], 1);

  //! @todo compute initial V2 = m - rX, SSE = yty - beta' * Xty
    
  glm_data.ErrorSS = yty - cblas_ddot(n, &beta[0], 1, &Xty[0], 1), 
  
  glm_data.V2 = geno_ind - rX;


  vector<double> SNPtSNP(geno_count), SNPty(geno_count);
  for (unsigned i=0; i<geno_count; ++i){
    //! these will never change for each SNP, so they could be moved out of all loops
    SNPtSNP[i] = cblas_ddot(geno_ind, &geno.values[i*geno_ind], 1, 
				&geno.values[i*geno_ind], 1);
    SNPty[i] = cblas_ddot(geno_ind, &geno.values[i*geno_ind], 1, &y[0], 1);
  }

  gettimeofday(&tstop, NULL);
  
  {
    const double computation_prep_time = 
      (tstop.tv_sec  - tstart.tv_sec) +
      (tstop.tv_usec - tstart.tv_usec)*1.e-6;
    cout << "Time required for computation prep: "
	 << computation_prep_time << " s" << endl;
  }
  
  gettimeofday(&tstart, NULL);

  // For each column of the geno array, set up the "X" matrix,
  // call the GLM routine, and store the computed p value.  Note
  // we are assuming y is a vector for now (because GLM currently expects
  // a vector for its second argument) but this could be generalized later.
  for (unsigned i=0; i<geno_count; ++i){

    // compute Xt * SNP    
    vector <double> XtSNP(geno_ind);

    //! @todo check m, n, lda
    //! @todo can be computed incrementally between major iterations
    cblas_dgemv(CblasColMajor, 
		CblasTrans,
		geno_ind,
		X.get_n_cols(),
		1.0,
		&X.values[0],
		geno_ind,
		&geno.values[i*geno_ind], // address of first SNP element
		1, // INCX = 1 because geno is column-major
		0.0,
		&XtSNP[0],
		1);

    // Call the glm function.  Note that X is currently overwritten by this function,
    // and therefore would need to be re-formed completely at each iteration...
    //glm(X, y, Kt, glm_data);
    /*
      ~200 us per SNP on Core i3, ~106 us per SNP on Core i7
     */
    //glm(X, XtX, XtXi, XtSNP, SNPtSNP, SNPty, yty, Kt, Xty, rX, glm_data);
    
    /*
      170 us per SNP on Core i3, 92 us per SNP on Core i7
     */
    plm(X, XtXi, XtSNP, SNPtSNP[i], SNPty[i], yty, Xty, rX, glm_data, 
	glm_data_new);

    // for p-val: p = 1 - fcdf(F, V1, V2), V1 = old V2 - new V2 (i.e. 0 or 1)
    // if V1 = 0, ignore; F is undefined
    // Store the computed value in an array
    Fval[i] = glm_data_new.F;
    V2s[i] = glm_data_new.V2; 
  }

  // Finish timing the computations
  gettimeofday(&tstop, NULL);

  {
    // Compute time taken for IO
    const double computation_elapsed_time = 
      (double)(tstop.tv_sec  - tstart.tv_sec) +
      (tstop.tv_usec - tstart.tv_usec)*1.e-6;
  
    cout << "Time required for computations: "
	 << computation_elapsed_time << " s" << endl;
    cout << "Time per SNP: " << computation_elapsed_time / geno_count 
	 << " s" << endl;
  }

  write("Fval.dat", Fval);
    
  return 0;
}


// On my Macbook (seconds) timings (algorithm assumes a Kt-matrix). 
// The algorithm which assumes a Kt-vector is only marginally faster, if at all...

// Avg. Time for I/O 
// (1.2631+1.24829+1.29038+1.26401+1.24963)/5 = 1.263082

// Avg. Computation time
// (0.279372+0.282701+0.271749+0.275394+0.274847)/5 = .2768126

// Percentage of total time used for computations:
// .2768126/(1.263082+.2768126) = .17976074
