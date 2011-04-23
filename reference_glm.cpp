/*! @todo
  - change input format to binary, or possibly zipped binary
  - get size of GPU RAM at run time
  - use GPU RAM size to estimate how many SNPs can be processed at a time
 */

// System header files
#include <sys/time.h>
#include <string>
#include <fstream>
#include <vector>
#include <limits>
#include <algorithm>
#include <stdlib.h>
#include <sstream>
#include <mpi.h>
#include <gsl/gsl_cdf.h>

// Local project includes
#include "fortran_matrix.h"
#include "glm.h"
#include "GetPot"
#include "svd.h"
#include "cblas.h"
#include "type.h"
#include "tvUtil.h"
#include "plm.h"

#define iterationLimit 50

const uint64_t readSize = 1024 * 1024 * 32;
const uint64_t readLength = readSize / sizeof(double);


using namespace std;

/*!
  divide total by ranks such that ranks * result >= total
 */
uint64_t adjust(uint64_t total, unsigned ranks){
  uint64_t result = total / ranks;
  if(total % ranks)
    result++;
  return result;
}

unsigned int nextPow2( unsigned int x ) {
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return ++x;
}

int readInputs(unsigned id, uint64_t myOffset, uint64_t mySize, string path, 
	       string fixed_filename, 
	       string geno_filename, 
	       string y_filename,
	       FortranMatrix &fixed, FortranMatrix &geno, vector<double> &y){
  // Read the "fixed" array from file.
  // assume this is small.
  {
    // Open "fixed" file for reading
    string filename = path + "/" + fixed_filename;
    ifstream fixed_file;
    fixed_file.open(filename.c_str(), ios::binary | ios::in);
    if (!fixed_file){
      cerr << "Failed to open fixed file: " << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    fixed_file.read((char*)&fixed.values[0], fixed.values.size() * sizeof(double));
    uint64_t temp = fixed.get_n_rows();
    fixed.n_rows = fixed.n_cols;
    fixed.n_cols = temp;
    fixed.transpose_self();
    for(unsigned row = 0; row < fixed.get_n_rows(); row++)
      for(unsigned col = 0; col < fixed.get_n_cols(); col++)
	if(fixed(row, col))
	  cout << row << " " << col << endl;
    fixed.writeD("fixed.dat");
  }    
  // Read the geno array from file
  /* each MPI rank reads a section of size
     total size / number of ranks
  */
  {

    // Open geno file for reading
    string filename = path + "/" + geno_filename;
    FILE *geno_file = fopen(filename.c_str(), "r");
    if(!geno_file){
      cerr << "failed to open geno file" << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if(fseeko(geno_file, myOffset, SEEK_SET)){
      cerr << "seek failed" << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    uint64_t readCount = 0, left = mySize / sizeof(double);
    while(left){
      size_t status = fread(&geno.values[readCount], sizeof(double), 
			    min(readLength, left), 
			    geno_file);
      if(status != min(readLength, left)){
	cerr << "read failed on id " << id << endl;
	if(feof(geno_file))
	  cerr << id << " eof" << endl;
	else if(ferror(geno_file))
	  cerr << id << " " << ferror(geno_file) << endl;
	MPI_Abort(MPI_COMM_WORLD, 1);
      }
      left -= status;
      readCount += status;
    }
    geno.writeD("geno.dat");
  }
  
  
  // Read the y-array from file.  Currently stored as a vector since that
  // is how it is passed to the glm function, but could be changed to a
  // FortranMatrix with one column...

  // assume this is small

  {
    // Open y file for reading
    ifstream y_file;
    y_file.open((path + "/" + y_filename).c_str(), ios::binary | ios::in);

    if (!y_file){
      cerr << "Failed to open file: " << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    y_file.read((char*)&y[0], y.size() * sizeof(double));
    writeD("y.dat", y);
  }
}

void compPrepare(FortranMatrix &X, FortranMatrix &fixed, unsigned fixed_count, FortranMatrix &XtX, vector<double> &Xty, vector<double> &y, FortranMatrix &U, vector<double> &S, FortranMatrix &Vt, unsigned &rX, vector<double> &beta, unsigned &n, double &tol, FortranMatrix &XtXi, double &yty, GLMData &glm_data, unsigned &geno_ind, unsigned &geno_count, const unsigned &m, FortranMatrix &geno, FortranMatrix &XtSNP, vector<double> &SNPty, vector<double> &SNPtSNP){
  // Fill first column of X with 1's
  for (unsigned i=0; i<X.get_n_rows(); ++i)
    X(i,0) = 1.;

  // Fill next fixed_count columns with the fixed array contents
  memcpy(&X.values[X.get_n_rows()], &fixed.values[0], 
	 fixed_count * X.get_n_rows() * sizeof(double));

  
  /*! precompute
    XtX
    XtXi
    rank(XtXi)
    
    for each SNP:
    SNPty
    SNPtSNP
   */
  
  XtX = matmat(X, X, true, false); 
  XtX.writeD("XtX.dat");

  // Solve (X^T * X)*beta = X^T*y for beta.  Note that X and X^T * X
  // have the same rank.

  // Initialize SVD components, A = U * S * V^T
  Xty = matvec(X, y, /*transX=*/true);
  
  writeD("Xty.dat", Xty);

  // Create the SVD of X^T * X 
  svd_create(XtX, U, S, Vt);
  U.writeD("U.dat");
  writeD("S.dat", S);

  Vt.writeD("Vt.dat");
  rX = svd_apply(U, S, Vt, /*result=*/beta, Xty);

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
  yty = cblas_ddot(y.size(), &y[0], 1, &y[0], 1);

  //! @todo compute initial V2 = m - rX, SSE = yty - beta' * Xty
    
  glm_data.ErrorSS = yty - cblas_ddot(n, &beta[0], 1, &Xty[0], 1), 
  
  glm_data.V2 = geno_ind - rX;

  
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


}

void compUpdate(vector<unsigned> &snpMask, unsigned &maxFIndex, FortranMatrix &X, 
		FortranMatrix &XtXi, FortranMatrix &XtSNP, 
		vector<double> &SNPtSNP, vector<double> &SNPty, double &yty,
		vector<double> &Xty, unsigned &rX, GLMData &glm_data,
		unsigned &n, unsigned &geno_count, const unsigned &m, 
		FortranMatrix &geno){
  snpMask[maxFIndex] = 1;

  glm(X, XtXi, &XtSNP(0, maxFIndex), SNPtSNP[maxFIndex], SNPty[maxFIndex], yty, Xty, 
      rX, glm_data);
  XtXi.writeD("XtXi.dat");
  XtSNP.resize_retain(n+1, geno_count);
  cblas_dgemv(CblasColMajor, CblasTrans, 
	      m, geno_count,
	      1.0, 
	      &geno.values[0],
	      m,
	      &geno(0, maxFIndex),
	      1,
	      0,
	      &XtSNP(n, 0),
	      n+1
	      );
  XtSNP.writeD("XtSNP.dat");


  // update host X
  X.resize_retain(m, n + 1);
  memcpy(&X(0, n), &geno(0, maxFIndex), m* sizeof(ftype));
  X.writeD("X.dat");
}

void getInputs(string &path, string &fixed_filename, string &geno_filename, 
	       string &y_filename, unsigned &fixed_count, unsigned &geno_ind,
	       unsigned &geno_count, int id){
  GetPot input_file("reference_glm.in");
  
  path = input_file("path", path.c_str());
  
  // File containing the "population structure".  It is a 4892-by-26 matrix
  fixed_filename = input_file("fixed_filename", "");
  if(fixed_filename == ""){
    if(!id)
      cerr << "invalid fixed_filename in input file" << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // File containing the genotypes.  It is a 4892-by-79 matrix.
  geno_filename = input_file("geno_filename", "");
  if(geno_filename == ""){
    if(!id)
      cerr << "invalid geno_filename in input file" << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // File containing the phenotypes.  It is a 4892-by-1 matrix.  The file is
  // arranged in a single column.
  y_filename = input_file("y_filename", "");
  if(y_filename == ""){
    if(!id)
      cerr << "invalid y_filename in input file" << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // In Matlab, these sizes are inferred from the data.  In C++, we hard-code them
  // to make reading the data simpler...
  fixed_count = input_file("fixed_count", 0); // rows, columns of the fixed array
  if(fixed_count == 0){
    if(!id)
      cerr << "invalid fixed_count in input file" << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  
  geno_ind = input_file("geno_ind", 0);
  geno_count = input_file("geno_count", 0); // rows, columns of the geno array
  if(geno_count == 0){
    if(!id)
      cerr << "invalid geno_count in input file" << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if(geno_ind == 0){
    if(!id)
      cerr << "invalid geno_ind in input file" << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

}

int main(int argc, char **argv)
{
  int id, numProcs;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  
  //! @todo this should be a command line parameter
  ftype entry_limit = 0.2;

  // Timing variables
  timeval tstart, tstop;

  //! @todo list of input files should be command line parameter

  // Create input file object.  Put the path to your data files here!
  // The path on my system to the location of the data files.  Don't forget the trailing
  // slash here, as this will be prepended to the filename below
  string path = "./"; // default path is the current directory.
  string fixed_filename;
  string geno_filename;
  string y_filename;
  unsigned fixed_count;
  unsigned geno_ind,
    geno_count; // rows, columns of the geno array
  
  getInputs(path, fixed_filename, geno_filename, y_filename, fixed_count, 
	    geno_ind, geno_count, id);

  const unsigned m = geno_ind;

  uint64_t totalSize = geno_count * geno_ind * sizeof(double);
  uint64_t perRankLength = adjust(geno_count, numProcs) * 
    geno_ind;
  uint64_t perRankSize =  perRankLength * sizeof(double);    
  uint64_t myOffset = id * perRankSize;
  uint64_t mySize = myOffset + perRankSize <= totalSize ?
    perRankSize : totalSize - myOffset;
  uint64_t mySNPs = mySize / sizeof(double) / geno_ind;

  // Matrix objects for storing the input data
  FortranMatrix fixed(geno_ind, fixed_count);
  FortranMatrix geno(geno_ind, mySNPs);
  vector<double> y(geno_ind);

  // Begin timing the file IO for all 3 files
  gettimeofday(&tstart, NULL);
  readInputs(id, myOffset, mySize, path, fixed_filename, geno_filename, y_filename, fixed, geno, y);
  gettimeofday(&tstop, NULL);
  
  cout << "Time required for I/O: " << tvDouble(tstop - tstart) << " s" << endl;
  
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
  FortranMatrix XtX;
  vector<double> S;
  FortranMatrix U, Vt, XtXi;

  GLMData glm_data, glm_data_new;
  vector<double> beta;
  vector<double> Xty;
  unsigned rX;
  double tol;
  double yty;
  vector<double> SNPtSNP(geno_count), SNPty(geno_count);;
  FortranMatrix XtSNP(n, geno_count);
  
  // Begin timing the computations
  gettimeofday(&tstart, NULL);

  compPrepare(X, fixed, fixed_count, XtX, Xty, y, U, S, Vt, rX, beta, n, tol, XtXi, 
	      yty, glm_data, geno_ind, geno_count, m, geno, XtSNP, SNPty, 
	      SNPtSNP);

  gettimeofday(&tstop, NULL);
  
  cout << "Time required for computation prep: "
       << tvDouble(tstop - tstart) << " s" << endl;
  
  ftype *d_snptsnp, *d_Xtsnp, *d_snpty,
    *d_f;
  //! @todo could use d_f also as a mask
  unsigned *d_snpMask;
  vector<unsigned> snpMask(geno_count, 0);
  // column-major with padding
  size_t d_XtsnpPitch;
  
  copyToDevice(geno_count, n, 
	       d_snptsnp, d_Xtsnp, d_XtsnpPitch, d_snpty, d_snpMask, d_f,
	       SNPtSNP, XtSNP, 
	       SNPty, Xty, XtXi, snpMask);
  
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
    
    unsigned maxFIndex;
    try{
      maxFIndex = plm_GPU(geno_count, n, n * sizeof(ftype), 
	m ,        
	d_snptsnp, 
	d_Xtsnp, 
	d_XtsnpPitch, 
	glm_data.ErrorSS, glm_data.V2, 
	d_snpty, 
	d_snpMask,
	d_f);
    } catch(int e){
      exit(e);
    }
    
    // for p-val: p = 1 - fcdf(F, V1, V2), V1 = old V2 - new V2 (i.e. 0 or 1)
    // if V1 = 0, ignore; F is undefined
    getMaxF(iteration, geno_count, Fval, maxFIndex, d_f);
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
  
    /*! update 
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
    compUpdate(snpMask, maxFIndex, X, XtXi, XtSNP, SNPtSNP, SNPty, yty, Xty, rX, glm_data, n , geno_count, m, geno);
    
    copyUpdateToDevice(geno_count, n, d_snpMask, maxFIndex, d_Xtsnp, 
		       d_XtsnpPitch, snpMask, XtSNP, XtXi, Xty);
  
    n++;

    {
      cout << "GPU time required for computations: "
	   << getGPUCompTime() << " s" << endl;
      cout << "GPU time per SNP: " << getGPUCompTime() / geno_count 
	   << " s" << endl;

      cout << "GPU time required for reduction: "
	   << getGPUMaxTime() << " s" << endl;
      cout << "GPU time per SNP: " << getGPUMaxTime() / geno_count 
	   << " s" << endl;
    }
    iteration++;
  } // while(1)

  MPI_Finalize();
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
