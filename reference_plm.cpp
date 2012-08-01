/*


Copyright (c) 2011, The Arizona Board of Regents on behalf of 
The University of Arizona

All rights reserved.

Developed by Peter Bailey, Tapasya Patki, and Greg Striemer with
support from the iPlant Collaborative as a collaboration between
participants at BIO5 at The University of Arizona (the primary hosting
institution), Cold Spring Harbor Laboratory, The University of Texas
at Austin, and individual contributors. Find out more at
http://www.iplantcollaborative.org/.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are
met:

 * Redistributions of source code must retain the above copyright 
   notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright 
   notice, this list of conditions and the following disclaimer in the 
   documentation and/or other materials provided with the distribution.
 * Neither the name of the iPlant Collaborative, BIO5, The University 
   of Arizona, Cold Spring Harbor Laboratory, The University of Texas at 
   Austin, nor the names of other contributors may be used to endorse or 
   promote products derived from this software without specific prior 
   written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


*/
/*! @todo print id numbers padded with 0s
  @todo use sched_setaffinity() to select a CPU core
  - if using two cores per node, use separate sockets

  @todo integrate SNP projection; 
  full SNPs are only used in a couple
  places (SNPtSNP and XtSNP prep, and XtSNP update).  
  If we were to reproject them when necessary, we could
  handle larger input sets.

  caveats: SNP projection currently operates per individual, not per
  SNP, producing row-major outputs.

  timing: Matlab on fermi does 3300 individuals per second with 4892
  SNPs (126081634 Bps), writing the results to disk.  Without writing
  the results to disk, the rate is 11000 (431572070 Bps). 

  I expect these number would be higher for a C implementation, as the
  input data is only 128 MB in ASCII format, uncompressed, and Matlab
  uses a single core.
 */

// System header files
#include <sys/time.h>
#include <string>
#include <string.h>
#include <fstream>
#include <vector>
#include <limits>
#include <algorithm>
#include <stdlib.h>
#include <sstream>
#include <iomanip>
#include <getopt.h>
#include <mpi.h>
#include <gsl/gsl_cdf.h>
#include <sys/sysinfo.h>

extern "C"{
#include <cblas.h>
}
// Local project includes
#include "fortran_matrix.h"
#include "glm.h"
#include "svd.h"

#include "type.h"
#include "tvUtil.h"
#include "plm.h"
#include "md5.h"

//! @todo match this to Lustre stripe size?
const uint64_t readSize = 1024 * 1024 * 32;
const uint64_t readLength = readSize / sizeof(double);
unsigned verbosity = 0;
unsigned iterationLimit = 50;

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

int readInputs(unsigned id, uint64_t myOffset, uint64_t mySize, 
	       uint64_t mySNPs,
	       string fixed_filename, 
	       string geno_filename, 
	       string y_filename,
	       FortranMatrix &fixed, FortranMatrix &geno, vector<double> &y,
	       bool rowMajor = false // applies only to geno for now
	       ){
  // Read the "fixed" array from file.
  // assume this is small.
  if(fixed_filename != ""){
    // Open "fixed" file for reading
    ifstream fixed_file;
    fixed_file.open(fixed_filename.c_str(), ios::binary | ios::in);
    if (!fixed_file){
      cerr << "Failed to open fixed file: " << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    //! assume row-major disk storage format for now
    fixed_file.read((char*)&fixed.values[0], fixed.values.size() * sizeof(double));
    fixed_file.close();
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
    FILE *geno_file = fopen(geno_filename.c_str(), "r");
    if(!geno_file){
      cerr << "failed to open geno file" << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if(fseeko(geno_file, myOffset, SEEK_SET)){
      cerr << "seek failed" << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    uint64_t readCount = 0, left = mySize / sizeof(double);
    /*! @todo convert from row-major format on disk to 
       column-major format in RAM
       - for now, just assume column-major on disk
     */
    if(rowMajor){
      double *array = (double*)malloc(mySNPs * sizeof(double));
      uint64_t row = 0;
      while(left){
#ifdef _DEBUG
	cout << "id " << id << " reading " << min(mySNPs, left) 
	     << " doubles from " << geno_filename << endl;
#endif
	size_t status = fread(array, sizeof(double), 
			      mySNPs, 
			      geno_file);
	if(status != mySNPs){
	  cerr << "read failed on id " << id << endl;
	  if(feof(geno_file))
	    cerr << id << " eof" << endl;
	  else if(ferror(geno_file))
	    cerr << id << " " << ferror(geno_file) << endl;
	  MPI_Abort(MPI_COMM_WORLD, 1);
	}
#ifdef _DEBUG
	for(int i = 0; i < mySNPs; i ++)
	  if(array[i] < 0){
	    cerr << "id " << id 
		 << " read a negative geno value at row " 
		 << row
		 << " column "
		 << myOffset / sizeof(double) + i
		 << "; aborting." 
		 << endl;
	    MPI_Abort(MPI_COMM_WORLD, 1);
	  }
#endif
	for(int i = 0; i < mySNPs; i ++)
	  geno(row, i) = array[i];
	left -= status;
	readCount += status;
	row++;
      }
      free(array);
    } else {
      while(left){
#ifdef _DEBUG
	cout << "id " << id << " reading " << min(readLength, left) 
	     << " doubles from " << geno_filename << endl;
#endif
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
#ifdef _DEBUG
	for(int i = readCount; i < readCount + status; i ++)
	  if(geno.values[i] < 0){
	    cerr << "id " << id 
		 << " read a negative geno value at global offset " 
		 << myOffset / sizeof(double) + i
		 << "; aborting." 
		 << endl;
	    MPI_Abort(MPI_COMM_WORLD, 1);
	  }
#endif
	left -= status;
	readCount += status;
      }
    }
#ifdef _DEBUG
    md5_state_t pms;
    md5_init(&pms);
    md5_append(&pms, (md5_byte_t*)&geno.values[0], mySize);
    md5_byte_t digest[16];
    md5_finish(&pms, digest);
    cout << "rank " << id << " md5 of geno:" ;
    for(int i = 0; i < 16; i++)
      cout << std::hex << (int)digest[i];
    cout << std::dec << endl;
      
#endif    
  }
  
  
  // Read the y-array from file.  Currently stored as a vector since that
  // is how it is passed to the glm function, but could be changed to a
  // FortranMatrix with one column...

  // assume this is small

  {
    // Open y file for reading
    ifstream y_file;
    y_file.open(y_filename.c_str(), ios::binary | ios::in);

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
		 FortranMatrix &geno, FortranMatrix &XtSNP, 
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

    //! @todo fix XtSNP lda for incremental computation
    /*! XtSNP is computed incrementally between major iterations.
      The initial computation requires more time, though, and
      @todo XtSNP could be computed as the geno data is read 
      from disk
     */
  cblas_dgemm(CblasColMajor,
	      CblasTrans,
	      CblasNoTrans,
	      n, 
	      mySNPs,
	      geno_ind,
	      1.0,
	      &X.values[0],
	      geno_ind, 
	      &geno.values[0],
	      geno_ind,
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
  /*! @todo SNPty could also be computed as the geno data is read 
    from disk
   */
  cblas_dgemv(CblasColMajor,
	      CblasTrans,
	      geno_ind,
	      mySNPs,
	      1.0,
	      &geno.values[0],
	      geno_ind,
	      &y[0],
	      1,
	      0.0,
	      &SNPty[0],
	      1
	      );

  /*! @todo SNPtSNP could be computed as the geno data is read from disk
   */
  for (uint64_t i=0; i<mySNPs; ++i){
    //! these will never change for each SNP
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

  timeval tGLMStart, tGLMStop, tResizeStop, tXtSNPStop;
  
  // update XtXi, Xty
  // output glm_data
  gettimeofday(&tGLMStart, 0);
  glm(id, iteration, n, XtXi, nextXtSNP, nextSNPtSNP, nextSNPty, yty, Xty,
	rX, glm_data);
  gettimeofday(&tGLMStop, 0);
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
  //! @todo this is probably slow
  XtSNP.resize_retain(n+1, mySNPs);

  gettimeofday(&tResizeStop, 0);

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

  gettimeofday(&tXtSNPStop, 0);
#ifdef _DEBUG
  {
    stringstream ss;
    ss << "XtSNP_" << iteration << "p_" << id << ".dat";
    XtSNP.writeD(ss.str());
  }
#endif
  if(verbosity > 1){
    /*
      10/24/2011, 1M 4892-SNPs
      GLM: 30us
      resize: 30ms
      XtSNP: 5s
     */
    cout << "id " << id 
	 << " GLM update time: " << tvDouble(tGLMStop - tGLMStart) << endl
	 << "id " << id 
	 << " resize time: " << tvDouble(tResizeStop - tGLMStop) << endl
	 << "id " << id
	 << " XtSNP time: " << tvDouble(tXtSNPStop - tResizeStop) << endl;
  }
}

template <class T> void write(const char *filename, const vector<T> &list){
  std::fstream file;
  file.open(filename, std::fstream::out);
  if(file.fail()){
    std::cerr << "failed to open file: " << filename << std::endl;
    exit(1);
  }
  file << std::setprecision(15);
  //file << std::scientific << std::setprecision(5);
  //  file << std::scientific << std::setprecision(15);
  
  unsigned m = list.size();

  for( int i = 0; i < m; i++ ) {
    file /*<< std::setw(13)*/ << " " << list[i] << endl;
  }
  file.close();
}


void printGlobalTime(timeval &tGlobalStart, timeval &tGlobalStop, 
		     double MPITime, double diskIOTime,
		     uint64_t mySNPs, unsigned iteration, unsigned id){
  gettimeofday(&tGlobalStop, NULL);
  
  cout << "id " << id << " total time: " 
       << tvDouble(tGlobalStop - tGlobalStart) + diskIOTime
       << "s" << endl;
  cout << "id " << id << " total time, not including disk/MPI IO: " 
       << tvDouble(tGlobalStop - tGlobalStart) - MPITime
       << "s" << endl;
  cout << "id " << id << " total time (-IO) per SNP per iteration: " 
       << (tvDouble(tGlobalStop - tGlobalStart) - MPITime) / mySNPs / iteration
       << "s" << endl;
}

void printUsage(char *name){
  cout << "usage: " << name << " [-f <fixed effects file> --num_fixed <number of fixed effects>] -g <SNP data file> --num_geno <number of SNPs> -r <residuals file> --num_r <number of residuals> [-c] [-v<verbosity level>] [-e SNP entry limit] [-l <max # of iterations>]" 
       << endl 
       << "where <input file> contains run-time settings" << endl;
}

int main(int argc, char **argv)
{
  int id, numProcs;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  ofstream logstream;
  streambuf *backupbuf, *logbuf;
  //redirect stdout
  stringstream ss_filename;
  ss_filename << "id_" << id << ".log";
  logstream.open(ss_filename.str().c_str());
  backupbuf = cout.rdbuf();
  logbuf = logstream.rdbuf();
  cout.rdbuf(logbuf);
  
  int opt;

  double entry_limit = 0.2;
  string fixed_filename = "";
  string geno_filename;
  string y_filename;
  unsigned fixed_count = 0,
    geno_ind; //rows 
  uint64_t geno_count; // columns of the geno array

  bool rowMajor = false;

  /*
    Skip a number of SNP entries for each individual; ideally, we will
    deal with the population and sample fields.
   */
  int skip = 0;

  int optIndex;
  struct option options[] = {{"num_fixed", required_argument, 0, 'n'},
			     {"num_geno", required_argument, 0, 'm'},
			     {"num_r", required_argument, 0, 'o'},
			     {"rowmajor", optional_argument, 0, 'p'},
			     {0,0,0,0}};
  
  /* get the following from the command line:

     -f <fixed effects filename>
     --num_fixed <number of fixed effects>

     -g <geno (SNP data) filename>
     --num_geno <number of SNPs>

     -r <residuals (y) filename>
     --num_r <number of residuals>

     -v <verbosity level>

     -e <entry limit>
  */

  while((opt = getopt_long(argc, argv, "cf:v::l:r:g:e:s:", options, &optIndex)) != -1){
    switch(opt){
    case 'e':
      entry_limit = atof(optarg);
      break;
    case 'c':
      break;
    case 'f':
      fixed_filename = optarg;
      break;
    case 'g':
      geno_filename = optarg;
      break;
    case 'r':
      y_filename = optarg;
      break;
    case 'v':
      if(optarg)
	verbosity = atoi(optarg);
      else
	verbosity = 1;
      break;
    case 'n':
      fixed_count = strtoul(optarg, 0, 0);
      break;
    case 'm':
      geno_count = strtoull(optarg, 0, 0);
      break;
    case 'o':
      geno_ind = strtoul(optarg, 0, 0);
      break;
    case 'l':
      iterationLimit = atoi(optarg);
      break;
    case 's':
      skip = strtoul(optarg, 0, 0);
      break;
    case 'p':
      rowMajor = true;
#ifdef _DEBUG
      cout << "Using row-major geno input" << endl;
#endif
      break;
    default:
      if(!id)
	printUsage(argv[0]);
      
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  if(geno_filename == "" || y_filename == ""){
    if(!id)
      printUsage(argv[0]);
    
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if((fixed_filename == "") ^ (!fixed_count)){
    if(!id)
      cout << "must supply both -f and --num_fixed, or neither" << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  timeval tstart, tstop;

  // global timing
  timeval tGlobalStart, tGlobalStop;

  //const unsigned m = geno_ind;

  uint64_t totalSize = geno_count * geno_ind * sizeof(double);
  uint64_t perRankSNPs = adjust(geno_count, numProcs);
  uint64_t perRankLength = perRankSNPs * geno_ind;
  uint64_t perRankSize =  perRankLength * sizeof(double);    
  uint64_t myOffset = id * perRankSize;
  uint64_t mySize = myOffset + perRankSize <= totalSize ?
    perRankSize : totalSize - myOffset;
  uint64_t mySNPs = mySize / sizeof(double) / geno_ind;
  uint64_t myStartSNP = id * perRankSNPs;

  if(verbosity > 1){
    cout << "id " << id << " has SNPs " << 
      myStartSNP << "-" << 
      myStartSNP + mySNPs - 1 << endl;
    
    struct sysinfo info;
    int status = sysinfo(&info);
    if(status){
      perror("error in sysinfo()");
    }
    
    cout << "id " << id << " SNPs require " 
	 << mySize / 1024.0 / 1024.0 / 1024.0 << " of " 
	 << info.totalram / 1024.0 / 1024.0 / 1024.0
	 << " GB host RAM" << endl;
    if(mySize > .95 * info.totalram){
      cerr << "id " << id << " SNP data is likely too large!" << endl;
    }
    if(!id)
      cout << "invoked as: ";
      for(int i = 0; i < argc; i++)
	cout << argv[i] << " ";
      cout << endl;
  }

  // Matrix objects for storing the input data
  FortranMatrix fixed(geno_ind, fixed_count);
  FortranMatrix geno(geno_ind, mySNPs);
  vector<double> y(geno_ind);
  vector<double> incomingSNP(geno_ind);
  vector<double> incomingXtSNP(iterationLimit + fixed_count + 1);
  double *nextSNP = &incomingSNP[0];
  double nextSNPtSNP, nextSNPty, *nextXtSNP = &incomingXtSNP[0];

  if(verbosity > 0)
    cout << "id " << id 
	 << " reading input data..." << endl;
  
  double 
    CPUCompUpdateTime = 0,
    MPITime = 0, diskIOTime = 0;

  // Begin timing the file IO for all 3 files
  gettimeofday(&tstart, NULL);
  readInputs(id, myOffset, mySize, mySNPs, fixed_filename, geno_filename, 
	     y_filename, fixed, geno, y, rowMajor);
  gettimeofday(&tstop, NULL);
  
  diskIOTime = tvDouble(tstop - tstart);
  if(verbosity > 0)
    cout << "id " << id 
	 << " I/O time: " << diskIOTime << " s" << endl;
  
  gettimeofday(&tGlobalStart, NULL);

  // Version B, Kt assumed a vector.
  //! assume Kt has a single non-zero element; a 1.0 at the end
  
  // only store global smallest P per iteration
  vector<double> Pval(iterationLimit); 

  // global indices of chosen SNPs from each iteration
  vector<int> chosenSNPs(iterationLimit);
  vector<int> chosenSNPsReduced(iterationLimit);

  // An array to hold the results of the GLM calculations
  vector<float> Fval(mySNPs);
 
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
	      yty, glm_data, geno_ind, mySNPs, geno, XtSNP, SNPty, 
	      SNPtSNP);

  gettimeofday(&tstop, NULL);
  
  if(verbosity > 0)
    cout << "id " << id 
	 << " computation prep time: "
	 << tvDouble(tstop - tstart) << " s" << endl;
  
  double *d_snptsnp, *d_Xtsnp, *d_snpty;
  float *d_f;
  //! @todo could use d_f also as a mask
  char *d_snpMask;
  vector<char> snpMask(mySNPs, 0);
  // column-major with padding
  size_t d_XtsnpPitch;
  
  // For each column of the geno array, set up the "X" matrix,
  // call the GLM routine, and store the computed p value.  Note
  // we are assuming y is a vector for now (because GLM currently expects
  // a vector for its second argument) but this could be generalized later.

    // Call the glm function.  Note that X is currently overwritten by this function,
    // and therefore would need to be re-formed completely at each iteration...
    /*
      ~200 us per SNP on Core i3, ~106 us per SNP on Core i7
     */
    
    /*
      170 us per SNP on Core i3, 92 us per SNP on Core i7
     */

  unsigned iteration = 0;
  while(iteration < iterationLimit && Pval[iteration] < entry_limit){
    // d_G, in constant memory
    // d_Xty, in constant memory
    
    /*! @todo this will need more bits when running more than 2^32 SNPs 
      on a single CPU process
    */
    int localMaxFIndex;

    float globalMaxF;

    // call CPU plm, get max F & index
    gettimeofday(&tstart, 0);
    int n = XtXi.get_n_rows();
      
    // G = XtXi

    // compute transpose of SNPtXG: 1xn
    vector<double> GtXtSNP(n, 0.0);
    for(uint64_t i = 0; i < mySNPs; i++){
      if(!snpMask[i]){
	int V2 = glm_data.V2;
	double ErrorSS = glm_data.ErrorSS;
	// previous V2 in glm_data

	//! @todo use cblas_dsymv for this
	cblas_dgemv(CblasColMajor,
		    CblasTrans, //! G is symmetric; but transpose is faster
		    n,
		    n,
		    1.0,
		    &XtXi.values[0],
		    n,
		    &XtSNP(0, i),
		    1,
		    0.0,
		    &GtXtSNP[0],
		    1);

	writeD("GtXtSNP.dat", GtXtSNP); // ok


	// compute SNPtXGXtSNP (scalar)
	// <SNPtX GtXtSNP> == <XtSNP GtXtSNP>
	double SNPtXGXtSNP = cblas_ddot(n, &XtSNP(0, i), 1, &GtXtSNP[0], 1);

	// compute S = Schur complement of partitioned matrix to invert
	double S = SNPtSNP[i] - SNPtXGXtSNP;
	if(S < doubleTol){ //! @todo if zero within tolerance
	  // bad news
	  Fval[i] = 0.0;
	  continue;
	}

	//S = 1.0 / S;

	// compute snpty - snptXGXty = snptMy == scalar
	// already know snpty, snptXG', Xty
	double SNPtMy = -cblas_ddot(n, &GtXtSNP[0], 1, &Xty[0], 1);
	SNPtMy += SNPty[i];

	double SSM = SNPtMy * SNPtMy / S;
	V2--;
	ErrorSS = ErrorSS - SSM;
	Fval[i] = V2 * SSM / ErrorSS;
      }else{
	Fval[i] = 0.0;
      }

      gettimeofday(&tstop, 0);
      double CPUCompTime = tvDouble(tstop - tstart);
      localMaxFIndex = cblas_isamax(mySNPs, &Fval[0], 1);
      gettimeofday(&tstart, 0);
      double CPUMaxTime = tvDouble(tstart - tstop);
      if(verbosity > 1){
	cout << "iteration " << iteration 
	     << " id " << id 
	     << " CPU computation time: "
	     << CPUCompTime << " s" << endl;
	cout << "iteration " << iteration 
	     << " id " << id 
	     << " CPU computation time per SNP: "
	     << CPUCompTime / mySNPs 
	     << " s" << endl;
	
	cout << "iteration " << iteration 
	     << " id " << id 
	     << " CPU reduction time: "
	     << CPUMaxTime << " s" << endl;
	cout << "iteration " << iteration 
	     << " id " << id 
	     << " CPU reduction time per SNP: "
	     << CPUMaxTime / mySNPs 
	     << " s" << endl;

      }
    }
    
    {
      stringstream ss;
      ss << "Fval_" << iteration << "_" << id << ".dat";
      writeD(ss.str(), Fval);
    }
    
    if(verbosity > 1){
      cout << "iteration " << iteration << " id " << id <<  
	" max F: " << Fval[localMaxFIndex] 
	   << " (local 0-index " << localMaxFIndex 
	   << ", global 0-index " << myStartSNP + localMaxFIndex << ")" << endl;
    }
    if(Fval[localMaxFIndex] <= 0){
      cerr << "error on iteration " << iteration << ": max F <= 0: " << Fval[localMaxFIndex] << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    gettimeofday(&tstart, NULL);
    // get max F value
    MPI_Allreduce(&Fval[localMaxFIndex], &globalMaxF, 1, MPI_FLOAT, MPI_MAX,
		  MPI_COMM_WORLD);
    gettimeofday(&tstop, NULL);
    MPITime += tvDouble(tstop - tstart);

    // get p value
    Pval[iteration] = 1 - gsl_cdf_fdist_P(globalMaxF, 1, glm_data.V2 - 1);

    if(Pval[iteration] > entry_limit){
      if(!id){
	cout << "p value (" << Pval[iteration] << ") > entry_limit (" << entry_limit 
	     << "); quitting" << endl;
	
	Pval.resize(iteration);
	chosenSNPs.resize(iteration);
      }
      break;
    }
      
    gettimeofday(&tstart, NULL);
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
    gettimeofday(&tstop, NULL);
    MPITime += tvDouble(tstop - tstart);

    if(verbosity > 1){
      if(!id)
	cout << "iteration " << iteration << " global max F on rank " << 
	  globalMinRankMaxF << ": " << globalMaxF << endl;
    }

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
      chosenSNPs[iteration] = localMaxFIndex + myStartSNP;
    }else{
      // receive SNP which yielded max F value
      nextSNP = &incomingSNP[0];
      nextXtSNP = &incomingXtSNP[0];
      localMaxFIndex = -1;
      chosenSNPs[iteration] = -1;
    }

    if(iteration + 1 >= iterationLimit){
      iteration++;
      break;
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
    gettimeofday(&tstart, NULL);
    MPI_Bcast(nextSNP, geno_ind, MPI_DOUBLE, globalMinRankMaxF,
	      MPI_COMM_WORLD);
    MPI_Bcast(nextXtSNP, n, MPI_DOUBLE, globalMinRankMaxF,
	      MPI_COMM_WORLD);
    MPI_Bcast(&nextSNPtSNP, 1, MPI_DOUBLE, globalMinRankMaxF,
	      MPI_COMM_WORLD);
    MPI_Bcast(&nextSNPty, 1, MPI_DOUBLE, globalMinRankMaxF,
	      MPI_COMM_WORLD);
    gettimeofday(&tstop, NULL);
    MPITime += tvDouble(tstop - tstart);
    
#ifdef _DEBUG
    {
      if(checkNeg(nextSNP, geno_ind))
	cerr << "rank " << id << " nextSNP negative at iteration " << iteration << endl;
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
      - Xty, (done)
      - geno, (done)
      - Xtsnp: matrix update from n x geno_count to n+1 x geno_count,
      - consists of appending (chosen SNP)' * geno = a new row
      - glm_data.ErrorSS, (done)
      - glm_data.V2, (done)
      - list of chosen SNP indices:
      To remove SNP from geno, set mask at SNP index.
     */

    gettimeofday(&tstart, NULL);
    compUpdate(id, iteration, XtXi, XtSNP, yty, Xty, rX, glm_data, n, mySNPs, geno_ind, 
	       geno, 
	       nextSNP, nextXtSNP, nextSNPtSNP, nextSNPty);
    gettimeofday(&tstop, NULL);
    CPUCompUpdateTime = tvDouble(tstop - tstart);

    if(verbosity > 1){
      cout << "iteration " << iteration 
	   << " id " << id 
	   << " CPU computation update time: "
	   << CPUCompUpdateTime << " s" << endl;
    }

    n++;
    iteration++;
  } // while(1)

  // reduce chosenSNPs to rank 0
  chosenSNPsReduced.resize(iteration);
  MPI_Reduce(&chosenSNPs[0], &chosenSNPsReduced[0], iteration, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

  if(iteration >= iterationLimit){
    if(!id)
      cout << "iteration limit (" << iterationLimit << ") reached" << endl;
  }

  if(verbosity > 0){
    printGlobalTime(tGlobalStart, tGlobalStop, MPITime, diskIOTime,
		    mySNPs, iteration, id);
    cout << "id " << id << " MPI communication time: " << MPITime << " s" 
	 << endl;
  }

  if(!id){
    write("Pval.dat", Pval);
    write("Pindices.dat", chosenSNPsReduced);
    cout << "iterations: " << iteration << endl;
  }

  // restore stdout
  cout.rdbuf(backupbuf);
  logstream.close();

  MPI_Finalize();
  return 0;
}
