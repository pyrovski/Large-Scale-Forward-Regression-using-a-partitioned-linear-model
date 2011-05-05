/*! @todo
  - get size of GPU RAM at run time
  - use GPU and host RAM size to estimate how many SNPs can be processed at a time
  -- host holds all SNPs plus derived data; currently, GPU holds only derived data
  - list of input files should be command line parameter
  - do comp update on GPU?
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
extern "C"{
#include <cblas.h>
}
// Local project includes
#include "fortran_matrix.h"
#include "glm.h"
#include "GetPot"
#include "svd.h"

#include "type.h"
#include "tvUtil.h"
#include "plm.h"

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
    
    //! assume row-major disk storage format for now
    fixed_file.read((char*)&fixed.values[0], fixed.values.size() * sizeof(double));
    fixed.transpose_dims();
    fixed.transpose_self();
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
    double *array = (double*)malloc(readSize);
    /*! @todo convert from row-major format on disk to 
       column-major format in RAM
       - for now, just assume column-major on disk
     */
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
#ifdef _DEBUG
    {
      stringstream ss;
      ss << "geno_" << id << ".dat";
      geno.writeD(ss.str());
    }
#endif
    free(array);
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
    if(!id)
      writeD("y.dat", y);
  }
  return 0;
}

void compPrepare(unsigned id, unsigned iteration, 
		 //FortranMatrix &X, 
		 FortranMatrix &fixed, 
		 unsigned fixed_count, FortranMatrix &XtX, vector<double> &Xty, 
		 vector<double> &y, FortranMatrix &U, vector<double> &S, 
		 FortranMatrix &Vt, unsigned &rX, vector<double> &beta, 
		 unsigned &n, double &tol, FortranMatrix &XtXi, double &yty, 
		 GLMData &glm_data, unsigned &geno_ind, uint64_t &mySNPs, 
		 const unsigned &m, FortranMatrix &geno, FortranMatrix &XtSNP, 
		 vector<double> &SNPty, vector<double> &SNPtSNP){

  FortranMatrix X(geno_ind/*4892*/, n/*27*/);

  // Fill first column of X with 1's


  for (unsigned i=0; i<X.get_n_rows(); ++i)
    X(i,0) = 1.;

  // Fill next fixed_count columns with the fixed array contents
  memcpy(&X(0, 1), &fixed.values[0], 
	 fixed_count * X.get_n_rows() * sizeof(double));

  
  /*! precompute
    XtX
    XtXi
    rank(XtXi)
    
    for each SNP:
    SNPty
    SNPtSNP
   */
  
  //XtX = matmat(X, X, true, false); 
  /*!
    construct XtX from X manually.
   */
  for(uint64_t col = 0; col < X.get_n_cols(); col++){
    if(col){
      cblas_dgemv(CblasColMajor,
		  CblasTrans,
		  geno_ind, 
		  col,
		  1.0,
		  &X.values[0],
		  geno_ind,
		  &X(0, col),
		  1,
		  0.0,
		  &XtX(0, col),
		  1);

      //! @todo if XtX is stored as symmetric matrix, we don't need this
      for(uint64_t row = 0; row < col; row++)
	XtX(col, row) = XtX(row, col);
    }
    XtX(col, col) = cblas_ddot(geno_ind, 
			       &X(0, col),
			       1,
			       &X(0, col),
			       1);
  }

  if(!id)
    XtX.writeD("XtX.dat");
  
  // Solve (X^T * X)*beta = X^T*y for beta.  Note that X and X^T * X
  // have the same rank.

  // Initialize SVD components, A = U * S * V^T
  Xty = matvec(X, y, /*transX=*/true);
  
  //! @todo this is altering XtX
  // Create the SVD of X^T * X 
  svd_create(XtX, U, S, Vt);

  rX = svd_apply(U, S, Vt, /*result=*/beta, Xty);

  if(!id){
    X.writeD("X.dat");
    XtX.writeD("XtX.dat");

    writeD("Xty_0.dat", Xty);

    U.writeD("U.dat");
    writeD("S.dat", S);
    Vt.writeD("Vt.dat");
    writeD("beta.dat", beta);
  }

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
#ifdef _DEBUG
  if(!id)
    {
      stringstream ss;
      ss << "XtXi_" << iteration << ".dat";
      XtXi.writeD(ss.str());
    }
#endif

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
	      mySNPs,
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
#ifdef _DEBUG
  {
    stringstream ss;
    ss << "XtSNP_" << iteration << "p_" << id << ".dat";
    XtSNP.writeD(ss.str());
  }
#endif

  //SNPty[i] = cblas_ddot(geno_ind, &geno.values[i*geno_ind], 1, &y[0], 1);
  cblas_dgemv(CblasColMajor,
	      CblasTrans,
	      m,
	      mySNPs,
	      1.0,
	      &geno.values[0],
	      m,
	      &y[0],
	      1,
	      0.0,
	      &SNPty[0],
	      1
	      );
  
  for (uint64_t i=0; i<mySNPs; ++i){
    //! these will never change for each SNP, so they could be moved out of all loops
    SNPtSNP[i] = cblas_ddot(geno_ind, &geno.values[i*geno_ind], 1, 
				&geno.values[i*geno_ind], 1);
  }
#ifdef _DEBUG
  {
    stringstream ss;
    ss << "SNPty_" << id << ".dat";
    writeD(ss.str(), SNPty);
  }
  {
    stringstream ss;
    ss << "SNPtSNP_" << id << ".dat";
    writeD(ss.str(), SNPtSNP);
  }
#endif


}

void compUpdate(unsigned id, unsigned iteration, 
		//FortranMatrix &X, 
		FortranMatrix &XtXi, FortranMatrix &XtSNP, 
		const double &yty,
		vector<double> &Xty, const unsigned &rX, GLMData &glm_data,
		const unsigned &n, const uint64_t &mySNPs, const unsigned &m, 
		FortranMatrix &geno,
		const double *nextSNP,
		const double *nextXtSNP,
		const double &nextSNPtSNP, 
		const double &nextSNPty){

  // update XtXi, Xty
  // output glm_data
  glm(id, iteration, n, XtXi, nextXtSNP, nextSNPtSNP, nextSNPty, yty, Xty, 
      rX, glm_data);
#ifdef _DEBUG
  if(!id)
    {
      {
	stringstream ss;
	ss << "XtXi_" << iteration << "p.dat";
	XtXi.writeD(ss.str());
      }
      {
	stringstream ss;
	ss << "Xty_" << iteration << "p.dat";
	writeD(ss.str(), Xty);
      }
    }
#endif
  XtSNP.resize_retain(n+1, mySNPs);
  cblas_dgemv(CblasColMajor, CblasTrans, 
	      m, mySNPs,
	      1.0, 
	      &geno.values[0],
	      m,
	      nextSNP,
	      1,
	      0,
	      &XtSNP(n, 0),
	      n+1
	      );
#ifdef _DEBUG
  {
    stringstream ss;
    ss << "XtSNP_" << iteration << "p_" << id << ".dat";
    XtSNP.writeD(ss.str());
  }
#endif

}

void getInputs(string input_filename, 
	       string &path, string &fixed_filename, string &geno_filename, 
	       string &y_filename, unsigned &fixed_count, unsigned &geno_ind,
	       uint64_t &geno_count, int id){

  // Create input file object.  Put the path to your data files here!
  GetPot input_file(input_filename.c_str());
  
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

  if(argc < 2){
    if(!id)
      cout << "usage: " << argv[0] << " <input file>" << endl 
	   << "where <input file> contains run-time settings" << endl;
    
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  
  //! @todo this should be a command line parameter
  double entry_limit = 0.2;

  // local timing
  timeval tstart, tstop;

  // global timing
  timeval tGlobalStart, tGlobalStop;

  // The path on my system to the location of the data files.  Don't forget the trailing
  // slash here, as this will be prepended to the filename below
  string path = "./"; // default path is the current directory.
  string fixed_filename;
  string geno_filename;
  string y_filename;
  unsigned fixed_count;
  unsigned geno_ind; //rows 
  uint64_t geno_count; // columns of the geno array
  string input_filename = argv[1];

  getInputs(input_filename, path, fixed_filename, geno_filename, y_filename, fixed_count, 
	    geno_ind, geno_count, id);

  const unsigned m = geno_ind;

  uint64_t totalSize = geno_count * geno_ind * sizeof(double);
  uint64_t perRankSNPs = adjust(geno_count, numProcs);
  uint64_t perRankLength = perRankSNPs * geno_ind;
  uint64_t perRankSize =  perRankLength * sizeof(double);    
  uint64_t myOffset = id * perRankSize;
  uint64_t mySize = myOffset + perRankSize <= totalSize ?
    perRankSize : totalSize - myOffset;
  uint64_t mySNPs = mySize / sizeof(double) / geno_ind;
  uint64_t myStartSNP = id * perRankSNPs;

#ifdef _DEBUG
  cout << "id " << id << " has SNPs " << 
    myStartSNP << "-" << 
    myStartSNP + mySNPs - 1 << endl;
#endif

  // Matrix objects for storing the input data
  FortranMatrix fixed(geno_ind, fixed_count);
  FortranMatrix geno(geno_ind, mySNPs);
  vector<double> y(geno_ind);
  vector<double> incomingSNP(geno_ind);
  vector<double> incomingXtSNP(iterationLimit + fixed_count + 1);
  double *nextSNP = &incomingSNP[0];
  double nextSNPtSNP, nextSNPty, *nextXtSNP = &incomingXtSNP[0];

  // Begin timing the file IO for all 3 files
  gettimeofday(&tGlobalStart, NULL);
  gettimeofday(&tstart, NULL);
  readInputs(id, myOffset, mySize, path, fixed_filename, geno_filename, 
	     y_filename, fixed, geno, y);
  gettimeofday(&tstop, NULL);
  
  cout << "id " << id 
       << " I/O time: " << tvDouble(tstop - tstart) << " s" << endl;
  
  // Version B, Kt assumed a vector.
  //! assume Kt has a single non-zero element; a 1.0 at the end
  
  // An array to hold the results of the GLM calculations
  double Pval;
  vector<double> Fval(mySNPs);
  vector<double> V2s(mySNPs);
  
  // Initialize the X-matrix.  The first column is all ones, the next
  // fixed_count columns are equal to the fixed matrix, and the last
  // column (which changes) is the i'th column of the geno array.
  unsigned n = fixed_count + 1;

  //! have to set size here; constructor for return value doesn't work in matmat
  FortranMatrix XtX(n, n);
  vector<double> S;
  FortranMatrix U, Vt, XtXi;

  GLMData glm_data, glm_data_new;
  vector<double> beta;
  vector<double> Xty;
  unsigned rX;
  double tol;
  double yty;
  vector<double> SNPtSNP(mySNPs), SNPty(mySNPs);
  FortranMatrix XtSNP(n, mySNPs);
  
  // Begin timing the computations
  gettimeofday(&tstart, NULL);

  compPrepare(id, 0, fixed, fixed_count, XtX, Xty, y, U, S, Vt, rX, 
	      beta, n, tol, XtXi, 
	      yty, glm_data, geno_ind, mySNPs, m, geno, XtSNP, SNPty, 
	      SNPtSNP);

  gettimeofday(&tstop, NULL);
  
  cout << "id " << id 
       << " computation prep time: "
       << tvDouble(tstop - tstart) << " s" << endl;
  
  double *d_snptsnp, *d_Xtsnp, *d_snpty,
    *d_f;
  //! @todo could use d_f also as a mask
  unsigned *d_snpMask;
  vector<unsigned> snpMask(mySNPs, 0);
  // column-major with padding
  size_t d_XtsnpPitch;
  
  double GPUCompTime = 0, GPUMaxTime = 0,
    GPUCopyTime = 0, GPUCopyUpdateTime = 0,
    CPUCompUpdateTime = 0;
  

  gettimeofday(&tstart, NULL);
  copyToDevice(mySNPs, n, 
	       d_snptsnp, d_Xtsnp, d_XtsnpPitch, d_snpty, d_snpMask, d_f,
	       SNPtSNP, XtSNP, 
	       SNPty, Xty, XtXi, snpMask);
  gettimeofday(&tstop, NULL);
  
  GPUCopyTime = tvDouble(tstop - tstart);
  
  cout << "id " << id 
       << " GPU copy time: "
       << GPUCopyTime << " s" << endl;
  cout << "GPU copy time per SNP: " << GPUCopyTime / mySNPs 
       << " s" << endl;

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
    
    if(iteration >= iterationLimit){
      if(!id)
	cout << "iteration limit (" << iterationLimit << ") reached" << endl;
      MPI_Finalize();
      return 0;
    }

    /*! @todo this will need more bits when running more than 2^32 SNPs 
      on a single CPU process
    */
    int localMaxFIndex;

    double globalMaxF;

    // ~3.5 us per SNP on Longhorn (FX5800)
    try{
      localMaxFIndex = plm_GPU(mySNPs, n, 
	m ,        
	d_snptsnp, 
	d_Xtsnp, 
	d_XtsnpPitch, 
	glm_data.ErrorSS, glm_data.V2, 
	d_snpty, 
	d_snpMask,
	d_f);
    } catch(int e){
      MPI_Abort(MPI_COMM_WORLD, e);
    }
    
    // for p-val: p = 1 - fcdf(F, V1, V2), V1 = old V2 - new V2 (i.e. 0 or 1)
    // if V1 = 0, ignore; F is undefined
    getMaxF(id, iteration, mySNPs, Fval, localMaxFIndex, d_f);
    cout << "iteration " << iteration << " id " << id <<  
      " max F: " << Fval[localMaxFIndex] 
	 << " (local 0-index " << localMaxFIndex 
	 << ", global 0-index " << myStartSNP + localMaxFIndex << ")" << endl;

    if(Fval[localMaxFIndex] <= 0){
      cerr << "error: max F <= 0: " << Fval[localMaxFIndex] << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // get max F value
    MPI_Allreduce(&Fval[localMaxFIndex], &globalMaxF, 1, MPI_DOUBLE, MPI_MAX,
		  MPI_COMM_WORLD);


    // get p value
    Pval = 1 - gsl_cdf_fdist_P(globalMaxF, 1, glm_data.V2 - 1);

    if(Pval > entry_limit){
      if(!id){
	cout << "p value (" << Pval << ") > entry_limit (" << entry_limit 
	     << "); quitting" << endl;
      }
      MPI_Finalize();
      return 0;
    }

    // determine a unique rank holding the max F value
    int globalMinRankMaxF;
    if(Fval[localMaxFIndex] == globalMaxF)
      MPI_Allreduce(&id, &globalMinRankMaxF, 1, MPI_INT, MPI_MIN, 
		    MPI_COMM_WORLD);
    else{
      int tempInt = numProcs + 1;
      MPI_Allreduce(&tempInt, &globalMinRankMaxF, 1, MPI_INT, MPI_MIN, 
		    MPI_COMM_WORLD);
    }

    if(!id)
      cout << "iteration " << iteration << " global max F on rank " << 
	globalMinRankMaxF << ": " << globalMaxF << endl;


    if(id == globalMinRankMaxF){
      // I have the max F value
      // send SNP which yielded max F value
      nextSNP = &geno(0, localMaxFIndex);
      nextXtSNP = &XtSNP(0, localMaxFIndex);
      nextSNPty = SNPty[localMaxFIndex];
      nextSNPtSNP = SNPtSNP[localMaxFIndex];
      snpMask[localMaxFIndex] = 1;
#ifdef _DEBUG
      cout << "iteration " << iteration << " id " << id 
	   << " masking index " << localMaxFIndex << endl;
#endif
    }else{
      // receive SNP which yielded max F value
      nextSNP = &incomingSNP[0];
      nextXtSNP = &incomingXtSNP[0];
      localMaxFIndex = -1;
    }

    /*
      need the following values from the SNP corresponding to the max F value:
      - SNP: vector(geno_ind)
      - SNPtSNP: scalar
      - SNPty: scalar
      - XtSNP: vector(n)
      
      for now, just broadcast them separately.
      it could be faster to send precalculated results in one transmission
      or recalculate them
     */
    MPI_Bcast(nextSNP, geno_ind, MPI_DOUBLE, globalMinRankMaxF,
	      MPI_COMM_WORLD);
    MPI_Bcast(nextXtSNP, n, MPI_DOUBLE, globalMinRankMaxF,
	      MPI_COMM_WORLD);
    MPI_Bcast(&nextSNPtSNP, 1, MPI_DOUBLE, globalMinRankMaxF,
	      MPI_COMM_WORLD);
    MPI_Bcast(&nextSNPty, 1, MPI_DOUBLE, globalMinRankMaxF,
	      MPI_COMM_WORLD);
    
#ifdef _DEBUG
    {
      stringstream ssSNP, ssXtSNP, ssSNPtSNP, ssSNPty;
      ssSNP << "nextSNP_" << id << "_" << iteration << ".dat";
      ssXtSNP << "nextXtSNP_" << id << "_" << iteration << ".dat";
      ssSNPtSNP << "nextSNPtSNP_" << id << "_" << iteration << ".dat";
      ssSNPty << "nextSNPty_" << id << "_" << iteration << ".dat";
      writeD(ssSNP.str(), nextSNP, geno_ind);
      writeD(ssXtSNP.str(), nextXtSNP, n);
      writeD(ssSNPtSNP.str(), &nextSNPtSNP, 1);
      writeD(ssSNPty.str(), &nextSNPty, 1);
    }
#endif

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

    gettimeofday(&tstart, NULL);
    compUpdate(id, iteration, XtXi, XtSNP, yty, Xty, rX, glm_data, n, mySNPs, m, 
	       geno, 
	       nextSNP, nextXtSNP, nextSNPtSNP, nextSNPty);
    gettimeofday(&tstop, NULL);
    CPUCompUpdateTime = tvDouble(tstop - tstart);

    gettimeofday(&tstart, NULL);
    copyUpdateToDevice(id, iteration, mySNPs, n, d_snpMask, localMaxFIndex, 
		       d_Xtsnp, d_XtsnpPitch, snpMask, XtSNP, XtXi, Xty);
    gettimeofday(&tstop, NULL);
    GPUCopyUpdateTime = tvDouble(tstop - tstart);
    
    n++;

    GPUCompTime = getGPUCompTime();
    GPUMaxTime = getGPUMaxTime();

    cout << "iteration " << iteration 
	 << " id " << id 
	 << " GPU computation time: "
	 << GPUCompTime << " s" << endl;
    cout << "iteration " << iteration 
	 << " id " << id 
	 << " GPU computation time per SNP: "
	 << GPUCompTime / mySNPs 
	 << " s" << endl;
    
    cout << "iteration " << iteration 
	 << " id " << id 
	 << " GPU reduction time: "
	 << GPUMaxTime << " s" << endl;
    cout << "iteration " << iteration 
	 << " id " << id 
	 << " GPU reduction time per SNP: "
	 << GPUMaxTime / mySNPs 
	 << " s" << endl;

    cout << "iteration " << iteration 
	 << " id " << id 
	 << " GPU copy update time: "
	 << GPUCopyUpdateTime << " s" << endl;

    cout << "iteration " << iteration 
	 << " id " << id 
	 << " CPU computation update time: "
	 << CPUCompUpdateTime << " s" << endl;

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
