// System header files
#include <sys/time.h>
#include <string>
#include <fstream>
#include <vector>
#include <limits>
#include <algorithm>
#include <cuda.h>
#include <cutil_inline.h>
#include <cublas.h>
#include <stdlib.h>
#include <sstream>


// Local project includes
#include "fortran_matrix.h"
#include "glm.h"
#include "GetPot"
#include "svd.h"
#include "cblas.h"

#define iterationLimit 50

#include "plm.cu"

using namespace std;

unsigned int nextPow2( unsigned int x ) {
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return ++x;
}

int main()
{
  //! @todo this should be a command line parameter
  ftype entry_limit = 0.2;

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
  unsigned y_ind = 4892;//, y_count = 1;        // rows, columns of the y vector

  const unsigned m = geno_ind;

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
	for (unsigned j=0; j<fixed.get_n_cols(); ++j)
	  fixed_file >> fixed(i,j);
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
	for (unsigned j=0; j<geno.get_n_cols(); ++j)
	  geno_file >> geno(i,j);
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
	for (unsigned i=0; i<y.size(); ++i)
	  y_file >> y[i];
    } else {
      cout << "Failed to open file!!" << endl;
      return 1;
    }
  }
  
  gettimeofday(&tstop, NULL);
  
  {
    // Compute time taken for IO
    const double io_elapsed_time = 
      (double)(tstop.tv_sec  - tstart.tv_sec) +
      (tstop.tv_usec - tstart.tv_usec)*1.e-6;
    
    cout << "Time required for I/O: " << io_elapsed_time << " s." << endl;
  }
  
  // Version A, Kt a general matrix.
  // Create the Kt matrix.  It has 1 row and fixed_count+2 columns.
  // The entries of Kt are all zero except for the last entry, which is 1.
  //FortranMatrix Kt(1,fixed_count+2);
  //Kt(0, fixed_count+1) = 1.; // Set last entry = 1

  // Version B, Kt assumed a vector.
  //! @todo assume Kt has a single non-zero element; a 1.0 at the end
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

  
  FortranMatrix XtSNP(n, geno_count);
    // compute Xt * SNP    

    //! @todo fix lda for incremental computation
    //! @todo can be computed incrementally between major iterations
  cblas_dgemm(CblasColMajor,
	      CblasTrans,
	      CblasNoTrans,
	      n, 
	      geno_count,
	      geno_ind,
	      1.0,
	      &X.values[0],
	      m, 
	      &geno.values[0],
	      m,
	      0.0,
	      &XtSNP.values[0],
	      n
	      );
  XtSNP.writeD("XtSNP.dat");

  vector<double> SNPtSNP(geno_count), SNPty(geno_count);
  //SNPty[i] = cblas_ddot(geno_ind, &geno.values[i*geno_ind], 1, &y[0], 1);
  cblas_dgemv(CblasColMajor,
	      CblasTrans,
	      m,
	      geno_count,
	      1.0,
	      &geno.values[0],
	      m,
	      &y[0],
	      1,
	      0.0,
	      &SNPty[0],
	      1
	      );
  writeD("SNPty.dat", SNPty);
  
  for (unsigned i=0; i<geno_count; ++i){
    //! these will never change for each SNP, so they could be moved out of all loops
    SNPtSNP[i] = cblas_ddot(geno_ind, &geno.values[i*geno_ind], 1, 
				&geno.values[i*geno_ind], 1);
  }

  gettimeofday(&tstop, NULL);
  
  {
    const double computation_prep_time = 
      (tstop.tv_sec  - tstart.tv_sec) +
      (tstop.tv_usec - tstart.tv_usec)*1.e-6;
    cout << "Time required for computation prep: "
	 << computation_prep_time << " s" << endl;
  }
  
  ftype *d_snptsnp, *d_Xtsnp, *d_snpty,
    *d_f;
  //! @todo could use d_f also as a mask
  unsigned *d_snpMask;
  vector<unsigned> snpMask(geno_count, 0);
  
  cutilSafeCall(cudaMalloc(&d_snpMask, geno_count * sizeof(unsigned)));
  cutilSafeCall(cudaMemcpy(d_snpMask, &snpMask[0], 
			   geno_count * sizeof(unsigned), 
			   cudaMemcpyHostToDevice));

  // column-major with padding
  size_t d_XtsnpPitch;
  
  //! @todo this won't be coalesced
  cutilSafeCall(cudaMalloc(&d_snptsnp, geno_count * sizeof(ftype)));
  cutilSafeCall(cudaMemcpy(d_snptsnp, &SNPtSNP[0], geno_count * sizeof(ftype), 
			   cudaMemcpyHostToDevice));

  cutilSafeCall(cudaMallocPitch(&d_Xtsnp, &d_XtsnpPitch, 
				(n + iterationLimit) * sizeof(ftype), 
				geno_count));
  cutilSafeCall(cudaMemcpy2D(d_Xtsnp, d_XtsnpPitch, &XtSNP.values[0], 
			     n * sizeof(ftype), n * sizeof(ftype), geno_count, 
			     cudaMemcpyHostToDevice));
  
  cutilSafeCall(cudaMalloc(&d_snpty, geno_count * sizeof(ftype)));
  cutilSafeCall(cudaMemcpy(d_snpty, &SNPty[0], geno_count * sizeof(ftype), 
			   cudaMemcpyHostToDevice));
  
  cutilSafeCall(cudaMemcpyToSymbol(d_G, &XtXi.values[0], n * n * sizeof(ftype)));
  cutilSafeCall(cudaMemcpyToSymbol(d_Xty, &Xty[0], n * sizeof(ftype)));
  
  cutilSafeCall(cudaMalloc(&d_f, geno_count * sizeof(ftype)));

  //gettimeofday(&tstart, NULL);
  cudaEvent_t start, stopKernel, stopMax;
  cudaEventCreate(&start);
  cudaEventCreate(&stopKernel);
  cudaEventCreate(&stopMax);
  
  
  // For each column of the geno array, set up the "X" matrix,
  // call the GLM routine, and store the computed p value.  Note
  // we are assuming y is a vector for now (because GLM currently expects
  // a vector for its second argument) but this could be generalized later.

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

    /*
    plm(X, XtXi, XtSNP, SNPtSNP[i], SNPty[i], yty, Xty, rX, glm_data, 
	glm_data_new);
    */
    
  unsigned iteration = 0;
  while(1){
    // d_G, in constant memory
    // d_Xty, in constant memory
    //unsigned threads = nextPow2(n);

    // clear last error
    cublasGetError();
    cudaEventRecord(start, 0);
    plm<<<geno_count, n, n * sizeof(ftype)>>>
      (m, 
       d_snptsnp, 
       d_Xtsnp, 
       d_XtsnpPitch, 
       glm_data.ErrorSS, glm_data.V2, 
       d_snpty, 
       d_snpMask,
       d_f);
    cudaEventRecord(stopKernel, 0);
    unsigned maxFIndex = cublasIdamax(geno_count, d_f, 1);
    cudaEventRecord(stopMax, 0);
    cutilSafeCall(cudaThreadSynchronize());
    // for p-val: p = 1 - fcdf(F, V1, V2), V1 = old V2 - new V2 (i.e. 0 or 1)
    // if V1 = 0, ignore; F is undefined
    // Store the computed value in an array
    //Fval[i] = glm_data_new.F;
    //V2s[i] = glm_data_new.V2; 
  
#ifndef _DEBUG
    cutilSafeCall(cudaMemcpy(&Fval[maxFIndex], &d_f[maxFIndex], sizeof(ftype),
			     cudaMemcpyDeviceToHost));
#else
    cutilSafeCall(cudaMemcpy(&Fval[0], d_f, geno_count * sizeof(ftype),
			     cudaMemcpyDeviceToHost));
    {
      stringstream ss;
      ss << "Fval_" << iteration << ".dat";
      writeD(ss.str(), Fval);
    }

#endif
    cout << "max F: " << Fval[maxFIndex] << " (" << maxFIndex << ")" << endl;
    
    // get p value
    if(Fval[maxFIndex] <= 0){
      // error
      cout << "error: max F <= 0" << endl;
      exit(1);
    }

    Pval[maxFIndex] = 1 - gsl_cdf_fdist_P(Fval[maxFIndex], 1, glm_data.V2 - 1);

    if(Pval[maxFIndex] > entry_limit){
      cout << "p value > entry_limit; quitting" << endl;
      break;
    }
  
    /*! @todo update 
      - X, (host only, done)
      - Xty, (done)
      - geno, (done)
      - Xtsnp: matrix update from n x geno_count to n+1 x geno_count,
      - consists of appending (chosen SNP)' * geno = a new row
      - glm_data.ErrorSS, (done)
      - glm_data.V2, (done)
      - list of chosen SNP indices:
      Assume data is in-core for GPU; i.e. don't recopy SNPs at each iteration.
      To remove SNP from geno, set mask at SNP index.
     */
    snpMask[maxFIndex] = 1;
    cutilSafeCall(cudaMemcpy(d_snpMask + maxFIndex, &snpMask[maxFIndex], 
			     sizeof(unsigned), cudaMemcpyHostToDevice));

    glm(X, XtXi, &XtSNP(0, maxFIndex), SNPtSNP[maxFIndex], SNPty[maxFIndex], yty, Xty, 
	rX, glm_data);
    XtSNP.resize_retain(n+1, geno_count);
    cblas_dgemv(CblasColMajor, CblasTrans, 
		n, geno_count, 
		1.0, 
		&XtSNP.values[0],
		n+1,
		&geno(0, maxFIndex),
		1,
		0,
		&XtSNP(n, 0),
		n+1
		);
    XtSNP.writeD("XtSNP.dat");

    //! copy updated XtSNP to GPU
    //! @todo not working?
    cutilSafeCall(
		  cudaMemcpy2D(d_Xtsnp + n, 
			       d_XtsnpPitch,
			       &XtSNP(n, 0), 
			       n + 1,
			       1,
			       geno_count,
			       cudaMemcpyHostToDevice)
		  );
    
    X.resize_retain(m, n + 1);
    memcpy(&X(0, n), &geno(0, maxFIndex), m);
    n++;

    {
      float computation_elapsed_time;
      cudaEventElapsedTime(&computation_elapsed_time, start, stopKernel);
      computation_elapsed_time /= 1000.0f;
      cout << "GPU time required for computations: "
	   << computation_elapsed_time << " s" << endl;
      cout << "GPU time per SNP: " << computation_elapsed_time / geno_count 
	   << " s" << endl;


      cudaEventElapsedTime(&computation_elapsed_time, stopKernel, stopMax);
      computation_elapsed_time /= 1000.0f;
      cout << "GPU time required for reduction: "
	   << computation_elapsed_time << " s" << endl;
      cout << "GPU time per SNP: " << computation_elapsed_time / geno_count 
	   << " s" << endl;
    }
    iteration++;
  } // while(1)
  
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
