#include <cuda.h>
#include <cutil_inline.h>
#include <cublas.h>
#include <sstream>
#include <iostream>

#include "type.h"
#include "fortran_matrix.h"
#include "reference_glm.h"

#include "plm.h"

#ifdef _DEBUG
#if __CUDA_ARCH__ >= 200
#define printGPU
#endif
#endif

__device__ int printBIDs(unsigned BID){
  return(BID == 0);
}

extern __shared__ double shared[];
extern __shared__ char sharedChar[];
#include "cuda_blas.cu"

using namespace std;
//__shared__ double fval; // scalar

__constant__ double d_Xty[fixedPlusIteration_limit + 1];

/*! extra space for padding.  Could make this triangular.
  If triangular, we could fit fixedPlusIteration_limit = 127, rather than 89
 */
__constant__ double d_G[(fixedPlusIteration_limit + 1)*(fixedPlusIteration_limit + 1)];

__global__ void plm(// inputs
		    const unsigned geno_count, // # of SNPs == # of blocks?
		    const unsigned m,          // rows of X
		    //const unsigned n,        // colums of X == threads/block
		    const double *snptsnp,      // scalar, unique to block
		    const double *Xtsnp,        // n x 1 vector, unique to block
		    const unsigned XtsnpPitchInWords, 
		    const double errorSS,       // scalar
		    const unsigned errorDF,       // scalar
		    //const double *G,          // symmetric n x n matrix in const mem
		    //const double *Xty,        // n x 1 vector in const mem
		    const double *snpty,        // scalar, unique to block
		    const char *snpMask,   // n x 1 vector
		    // outputs
		    float *f){
  /*! @todo could compute two SNPs per thread block.  
    This would ease the limitation of 8 thread blocks/MP for SM 1.3 devices,
    but might need some thread padding for warps.
   */

  double GtXtsnp; // each thread stores one element of each array // Xtsnp
  //! @todo these might use fewer registers if kept in shared memory
  double snptmy; // scalar

  double s; // scalar

  unsigned BID = blockIdx.x + gridDim.x * blockIdx.y;
  unsigned TID = threadIdx.x;

  #ifdef printGPU
  if(printBIDs(BID) && !TID)
    printf("BID: %u\n", BID);
  #endif

  if(BID >= geno_count)
    return;
  if(!TID)
    sharedChar[0] = snpMask[BID];

  __syncthreads();
  if(sharedChar[0]){
    // don't compute new F
    if(!TID)
      f[BID] = 0;
    return; 
  }
  // snptsnp - snptXGXtsnp:

  // 3n ops
  // 1 read per thread = n total
  double myXtsnp = *(Xtsnp + BID * XtsnpPitchInWords + TID);
  
  // GtXtsnp: n * (1 fmad, 1 imad?) per thread, 2n^2 total operations,
  // n^2 total reads
  GtXtsnp = vecRMatCSq(TID, BID, myXtsnp, blockDim.x, d_G, 
		       blockDim.x);  //! length of column plus padding (no padding)
  
  // snptsnp - snptXGXtsnp
  // ~2n operations, no reads
  dotRR(TID, blockDim.x, GtXtsnp, myXtsnp);

  // n operations
  s = snptsnp[BID] - shared[0];
#ifdef printGPU
  if(printBIDs(BID)){
    for(int i = 0; i < blockDim.x; i++){
      if(i == TID){
	printf("b%03u\tt%03u\tXtsnp: %1.10le\n", BID, TID, Xtsnp[BID * XtsnpPitchInWords + TID]);
	for(int j = 0; j < blockDim.x; j++)
	  printf("b%03u\tt%03u\tG[%d,%d]: %1.10le\n", BID, TID, i, j, d_G[i*blockDim.x + j]);
	printf("b%03u\tt%03u\tGtXtsnp: %1.10le\n", BID, TID, GtXtsnp);
	if(!TID){
	  printf("b%03u\tt%03u\tsnptsnp: %1.10le\n", BID, TID, snptsnp[BID]);
	  printf("b%03u\tt%03u\tsnptXGXtsnp: %1.10le\n", BID, TID, shared[0]);
	  printf("b%03u\tt%03u\ts: %1.10le\n", BID, TID, s);
	}
      }
      __syncthreads();
    }
  }
#endif
  // 1/(above)
  // n ops
  if(s > doubleTol){
    // snptmy
    // ~2n ops, n reads
    dotRG(TID, blockDim.x, GtXtsnp, d_Xty);
    
    if(!TID){
      // 1 op
      //s = (double)1/s;
      
      // 1 op
      //snptmy = -shared[0];
#ifdef printGPU
      if(printBIDs(BID)){
	printf("b%03u\tt%03u\tsnptXGXty: %1.10le\n", BID, TID, -shared[0]);
	printf("b%03u\tt%03u\tsnpty: %1.10le\n", BID, TID, snpty[BID]);
      }
#endif
      
      // 1 op
      snptmy = snpty[BID] - shared[0];

      // 2 ops
      double modelSS = (snptmy * snptmy) / s;

      // 1 op
      double errorSS2 = errorSS - modelSS;

      // 1 op
      unsigned V2 = errorDF - 1;

      // 2 ops
      f[BID] = modelSS / errorSS2 * V2;
#ifdef printGPU
      if(f[BID] < 0){
	printf("b%03u\tt%03u\tnew f: %1.10le < 0!\n", BID, TID, f[BID]);
	printf("b%03u\tt%03u\tnew modelSS: %1.10le\n", BID, TID, modelSS);
	printf("b%03u\tt%03u\terrorSS: %1.10le\n", BID, TID, errorSS2);
	printf("b%03u\tt%03u\tnew errorSS: %1.10le\n", BID, TID, errorSS2);
	printf("b%03u\tt%03u\tnew V2: %u\n", BID, TID, V2);
      }
      if(printBIDs(BID)){
	printf("b%03u\tt%03u\tmodelSS: %1.10le\n", BID, TID, modelSS);
	printf("b%03u\tt%03u\tnew errorSS: %1.10le\n", BID, TID, errorSS2);
	printf("b%03u\tt%03u\tnew V2: %u\n", BID, TID, V2);
	printf("b%03u\tt%03u\tf: %1.10le\n", BID, TID, f[BID]);
      }
#endif
      
    }
    return;
  } else {
    if(!TID)
      f[BID] = 0;
    return;
  }

  /*
    
    Selective GPU:
    ops: $3n+2n^2+2n+n+2n+1+2+1+1+2  =  2n^2+8n+7$
    
    reads: $n+n^2+n  =  n^2+2n, n^2$ of which are from constant memory
    
    writes: 1
  */
}

cudaEvent_t start, stopKernel, stopMax;

void initGrid(dim3 &grid, unsigned geno_count) throw(int){
  static unsigned old_geno_count = 0;
  static dim3 oldGrid;
  if(old_geno_count == geno_count){
    grid = oldGrid;
    return;
  }
  grid.x = geno_count;
  grid.y = 1;
  grid.z = 1;
  while(grid.x > 65535){
    grid.x += grid.x % 2;
    grid.x /= 2;
    grid.y *= 2;
  }
  // could probably do a better job factoring here; instead bail
  if(grid.y > 65535)
    throw(2);
  oldGrid = grid;
#ifdef _DEBUG
  cout << "grid: " << grid.x << "x" << grid.y << endl;
#endif
}

unsigned plm_GPU(unsigned geno_count, unsigned blockSize, 
		 unsigned m, double* d_snptsnp, double* d_Xtsnp, 
		 unsigned d_XtsnpPitch, double ErrorSS, unsigned V2, 
		 double* d_snpty, char* d_snpMask, float* d_f,
		 vector<float> &Fval) throw(int)
{
    cublasGetError();
#ifdef _DEBUG
    cutilSafeCall(cudaThreadSynchronize());
#endif
    cudaEventRecord(start, 0);
    dim3 grid;
    initGrid(grid, geno_count);
    
    plm<<<grid, blockSize, blockSize * sizeof(double)>>>
      (geno_count,
       m,
       d_snptsnp, 
       d_Xtsnp, 
       d_XtsnpPitch / sizeof(double), 
       ErrorSS, V2, 
       d_snpty, 
       d_snpMask,
       d_f);
    cudaEventRecord(stopKernel, 0);

#ifdef _DEBUG
    cutilSafeCall(cudaThreadSynchronize());
    cutilSafeCall(cudaMemcpy(&Fval[0], d_f, geno_count * sizeof(float),
			     cudaMemcpyDeviceToHost));
#endif

    cublasStatus status = cublasGetError();

    // cublas uses 1-based index
    int maxFIndex = cublasIsamax(geno_count, d_f, 1);
    cudaEventRecord(stopMax, 0);
    status = cublasGetError();
    if(status != CUBLAS_STATUS_SUCCESS){
      cerr << "cublas error in cublasIdamax(): " << status << endl;
      throw(1);
    }
    if(maxFIndex <= 0){
      cerr << "maxFIndex <= 0:" << maxFIndex << endl;
      throw(1);
    }
    maxFIndex -= 1;

    return maxFIndex;
}

/*!
  should only be called once
 */
int copyToDevice(const unsigned id, 
		 const unsigned verbosity,
		 const unsigned geno_count, const unsigned n, 
		 double *&d_snptsnp, double *&d_Xtsnp, size_t &d_XtsnpPitch, 
		 double *&d_snpty, char *&d_snpMask, float *&d_f,
		 const vector<double> &SNPtSNP, const FortranMatrix &XtSNP,
		 const vector<double> &SNPty,
		 const vector<double> &Xty, const FortranMatrix &XtXi, 
		 const vector<char> &snpMask){

  uint64_t snpMaskSize = geno_count * sizeof(char), 
    snptsnpSize = geno_count * sizeof(double), 
    XtsnpSize, 
    snptySize = geno_count * sizeof(double), 
    fSize = geno_count * sizeof(float);

  uint64_t totalSize;

  cudaError_t cudaStatus;
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  int device = deviceCount > 1 ? (id % 2) : 0;
  /*! @todo set device numbers from sge wayness parameter;
    for now, assume that if ID is even, use first GPU, 
    otherwise, use second GPU
  */
  cudaStatus = cudaSetDevice(device);
  if(cudaStatus != cudaSuccess){
    cerr << "id " << id << " error in cudaSetDevice()" << endl;
    return -1;
  }
  if(verbosity > 1){
    cout << "id " << id << " using GPU " << device + 1 << " of " << deviceCount
	 << endl;
  }
  
  cutilSafeCall(cudaMallocPitch(&d_Xtsnp, &d_XtsnpPitch, 
				(n + iterationLimit) * sizeof(double), 
				geno_count));
  XtsnpSize = d_XtsnpPitch * geno_count;
  totalSize = fSize + XtsnpSize + snptsnpSize + snpMaskSize + snptySize;

  cudaFree(&d_Xtsnp);

  struct cudaDeviceProp prop;
  cudaStatus = cudaGetDeviceProperties(&prop, device);
  if(cudaStatus != cudaSuccess){
    cerr << "id " << id << " error in cudaGetDeviceProperties()" << endl;
    return -1;
  }
  if(verbosity > 1)
    cout << "id " << id << " requires " << totalSize << "/" 
	 << prop.totalGlobalMem << " bytes global memory"  
	 << " (" << (100.0 * totalSize) / prop.totalGlobalMem 
	 << "%)"
	 << endl;
  if(totalSize >= prop.totalGlobalMem){
    cerr << "id " << id << " insufficient device memory" << endl;
    return -1;
  }

  /* this is checked by the compiler
  if(totalConstantSize >= prop.totalConstMem){
    cerr << "id " << id << " insufficient device constant memory" << endl;
    return -1;
  }
  */
  
  cutilSafeCall(cudaMallocPitch(&d_Xtsnp, &d_XtsnpPitch, 
				(n + iterationLimit) * sizeof(double), 
				geno_count));

  cutilSafeCall(cudaMalloc(&d_snpMask, geno_count * sizeof(char)));
  cutilSafeCall(cudaMemcpy(d_snpMask, &snpMask[0], 
			   geno_count * sizeof(char), 
			   cudaMemcpyHostToDevice));
  
  //! @todo this won't be coalesced
  cutilSafeCall(cudaMalloc(&d_snptsnp, geno_count * sizeof(double)));
  cutilSafeCall(cudaMemcpy(d_snptsnp, &SNPtSNP[0], geno_count * sizeof(double), 
			   cudaMemcpyHostToDevice));

  cutilSafeCall(cudaMemcpy2D(d_Xtsnp, d_XtsnpPitch, &XtSNP.values[0], 
			     n * sizeof(double), n * sizeof(double), geno_count, 
			     cudaMemcpyHostToDevice));
  
  cutilSafeCall(cudaMalloc(&d_snpty, geno_count * sizeof(double)));
  cutilSafeCall(cudaMemcpy(d_snpty, &SNPty[0], geno_count * sizeof(double), 
			   cudaMemcpyHostToDevice));
  
  cutilSafeCall(cudaMemcpyToSymbol(d_G, &XtXi.values[0], n * n * sizeof(double)));
  cutilSafeCall(cudaMemcpyToSymbol(d_Xty, &Xty[0], n * sizeof(double)));
  
  cutilSafeCall(cudaMalloc(&d_f, geno_count * sizeof(float)));

  cublasStatus status = cublasInit();
  if(status != CUBLAS_STATUS_SUCCESS){
    cerr << "id " << id << " error in cublasInit()" << endl;
    return -1;
  }
  cudaEventCreate(&start);
  cudaEventCreate(&stopKernel);
  cudaEventCreate(&stopMax);
  return 0;
}

void copyUpdateToDevice(unsigned id, unsigned iteration,  
			unsigned geno_count, unsigned n,
		       char *d_snpMask, 
		       int maxFIndex, double *d_Xtsnp, 
		       size_t d_XtsnpPitch,
		       const vector<char> &snpMask,
		       FortranMatrix &XtSNP, const FortranMatrix &XtXi,
		       const vector<double> &Xty){

  if(maxFIndex >= 0){
#ifdef _DEBUG
      cout << "iteration " << iteration << " id " << id 
	   << " masking index " << maxFIndex << " on device" 
	   << endl;
#endif
    cutilSafeCall(cudaMemcpy(d_snpMask + maxFIndex, &snpMask[maxFIndex], 
			     sizeof(char), cudaMemcpyHostToDevice));
#ifdef _DEBUG
    unsigned maskVal;
    cutilSafeCall(cudaMemcpy(&maskVal, d_snpMask + maxFIndex,
			     sizeof(char), cudaMemcpyDeviceToHost));
      cout << "iteration " << iteration << " id " << id 
	   << " mask index " << maxFIndex << ": "
	   << maskVal << endl;
#endif
  }

    //! copy updated XtSNP to GPU
    cutilSafeCall(cudaMemcpy2D(d_Xtsnp + n, 
			       d_XtsnpPitch,
			       &XtSNP(n, (unsigned)0), 
			       (n + 1) * sizeof(double),
			       sizeof(double),
			       geno_count,
			       cudaMemcpyHostToDevice));

    // update GPU G (const mem)
    cutilSafeCall(cudaMemcpyToSymbol(d_G, &XtXi.values[0], 
				     (n + 1) * (n + 1) * sizeof(double)));

    // update GPU Xty (const mem)
    // call fails unless we update the whole thing
    cutilSafeCall(cudaMemcpyToSymbol(d_Xty, &Xty[0], (n + 1) * sizeof(double)));
}

float getGPUCompTime(){
  float computation_elapsed_time;
  cudaEventElapsedTime(&computation_elapsed_time, start, stopKernel);
  return computation_elapsed_time / 1000.0f;
}

float getGPUMaxTime(){
  float computation_elapsed_time;
  cudaEventElapsedTime(&computation_elapsed_time, stopKernel, stopMax);
  return computation_elapsed_time / 1000.0f;
}

void getMaxFGPU(unsigned id, unsigned iteration, unsigned geno_count, 
	     vector<float> &Fval, 
	     unsigned maxFIndex, float *d_f){
#ifndef _DEBUG
    cutilSafeCall(cudaMemcpy(&Fval[maxFIndex], &d_f[maxFIndex], sizeof(float),
			     cudaMemcpyDeviceToHost));
#endif
}
