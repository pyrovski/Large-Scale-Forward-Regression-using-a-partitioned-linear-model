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
  std::cout << "grid: " << grid.x << "x" << grid.y << std::endl;
#endif
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
//! @todo should y be in constant mem?  Not if we use cublas for dgemv (SNPty)

__global__ void plm(// inputs
		    const unsigned geno_count, // # of SNPs == # of blocks?
		    const unsigned m,          // rows of X
		    //const unsigned n,        // colums of X == threads/block
		    const double *snptsnp,      // scalar, unique to block
		    const double *Xtsnp,        // n x 1 vector, unique to block
		    const unsigned XtsnpPitchInWords, 
		    const double errorSS,       // scalar
		    const unsigned errorDF,       // scalar
		    //const double *G,          // symmetric matrix in const mem
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
  if(!TID && printBIDs(BID))
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

  //#error fixme double read
  double myXtsnp = *(Xtsnp + BID * XtsnpPitchInWords + TID);
  // GtXtsnp
  GtXtsnp = vecRMatCSq(TID, myXtsnp, blockDim.x, d_G, 
		       blockDim.x);  //! length of column plus padding (no padding) 
  
  // snptsnp - snptXGXtsnp
  dotRR(TID, blockDim.x, GtXtsnp, myXtsnp);
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
  if(s > doubleTol){
    s = (double)1/s;
    
    // snptmy
    dotRG(TID, blockDim.x, GtXtsnp, d_Xty);
    snptmy = -shared[0];
    
    if(!TID){
#ifdef printGPU
      if(printBIDs(BID)){
	for(int i = 0; i < blockDim.x; i++){
	  if(i == TID){
	    printf("b%03u\tt%03u\tsnptXGXty: %1.10le\n", BID, TID, snptmy);
	    printf("b%03u\tt%03u\tsnpty: %1.10le\n", BID, TID, snpty[BID]);
	  }
	  __syncthreads();
	}
      }
#endif
      snptmy += snpty[BID];
      //! @todo this could be single precision
      float modelSS = snptmy * snptmy * s;

      double errorSS2 = errorSS - (double)modelSS;

      unsigned V2 = errorDF - 1;
      f[BID] = modelSS / errorSS2 * V2;
#ifdef printGPU
  if(printBIDs(BID)){
    for(int i = 0; i < blockDim.x; i++){
      if(i == TID){
	
	printf("b%03u\tt%03u\tmodelSS: %1.10le\n", BID, TID, modelSS);
	printf("b%03u\tt%03u\tnew errorSS: %1.10le\n", BID, TID, errorSS2);
	printf("b%03u\tt%03u\tnew V2: %u\n", BID, TID, V2);
	printf("b%03u\tt%03u\tf: %1.10le\n", BID, TID, f[BID]);
      }
      __syncthreads();
    }
  }
#endif

    }
    return;
  } else {
    if(!TID)
      f[BID] = 0;
    return;
  }
}

cudaEvent_t start, stopKernel, stopMax;

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

//! @todo convert for gpu
/*!
  should only be called once
 */
int copyToDevice(const unsigned id, // MPI rank
		 const unsigned verbosity,
		 const unsigned geno_count, // # of SNPs
		 const unsigned n, // # of columns in X
		 const unsigned geno_ind, // # length of each SNP
		 double *&d_snptsnp, double *&d_Xtsnp, size_t &d_XtsnpPitch, 
		 double *&d_snpty, char *&d_snpMask, float *&d_f, 
		 double *&d_geno, size_t &d_genoPitch, 
		 double *&d_Xt, size_t &d_XtPitch,
		 double *&d_y,
		 const vector<double> &Xty, 
		 const FortranMatrix &XtXi, // G
		 const FortranMatrix &Xt, // only for XtSNP
		 const FortranMatrix &geno,
		 double *&d_nextSNP,
		 const double *y){

  uint64_t snpMaskSize = geno_count * sizeof(char), 
    snptsnpSize = geno_count * sizeof(double), 
    XtsnpSize, 
    snptySize = geno_count * sizeof(double), 
    fSize = geno_count * sizeof(float),
    genoSize, // new
    XtSize, // new
    nextSNPSize = geno_ind * sizeof(double); // new

  uint64_t totalSize;

  
  cutilSafeCall(cudaMallocPitch(&d_Xtsnp, &d_XtsnpPitch, 
				(n + iterationLimit) * sizeof(double), 
				geno_count));
  XtsnpSize = d_XtsnpPitch * geno_count;
   
  cutilSafeCall(cudaMallocPitch(&d_geno, &d_genoPitch, geno_ind * sizeof(double), geno_count));
  if(d_genoPitch < geno_ind){
    cerr << "error allocating d_geno" << endl;
    return -1;
  }
  genoSize = geno_count * d_genoPitch;

  //! store Xt
  cutilSafeCall(cudaMallocPitch(&d_Xt, &d_XtPitch, n, geno_ind));
  if(d_XtPitch < n){
    cerr << "error allocating d_Xt" << endl;
    return -1;
  }
  XtSize = geno_ind * d_XtPitch;

  // Xty and G are in constant memory; compiler checks this
  totalSize = fSize + XtsnpSize + snptsnpSize + snpMaskSize + snptySize
    + genoSize + XtSize + nextSNPSize;

  struct cudaDeviceProp prop;
  cudaError_t cudaStatus;
  int device;
  cudaGetDevice(&device);
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

  cutilSafeCall(cudaMemcpy2D(d_geno, d_genoPitch, &geno.values[0], 
			     geno_ind * sizeof(double), 
			     geno_ind * sizeof(double),
			     geno_count, cudaMemcpyHostToDevice));
  cutilSafeCall(cudaMemcpy2D(d_Xt, d_XtPitch, &Xt.values[0],
			     n * sizeof(double), 
			     n * sizeof(double),
			     geno_ind, cudaMemcpyHostToDevice));
  cutilSafeCall(cudaMalloc(&d_snpMask, geno_count * sizeof(char)));

  cutilSafeCall(cudaMemsetAsync(d_snpMask, 0, geno_count * sizeof(char)));
  
  //! @todo these won't be coalesced; each block reads one value
  cutilSafeCall(cudaMalloc(&d_snptsnp, geno_count * sizeof(double)));  
  cutilSafeCall(cudaMalloc(&d_snpty, geno_count * sizeof(double)));
  
  cutilSafeCall(cudaMemcpyToSymbol(d_G, &XtXi.values[0], 
				   n * n * sizeof(double)));
  cutilSafeCall(cudaMemcpyToSymbol(d_Xty, &Xty[0], n * sizeof(double)));
  
  cutilSafeCall(cudaMalloc(&d_f, geno_count * sizeof(float)));
  cutilSafeCall(cudaMalloc(&d_nextSNP, geno_ind * sizeof(double)));

  cutilSafeCall(cudaMalloc(&d_y, geno_ind * sizeof(double)));
  cutilSafeCall(cudaMemcpy(d_y, y, geno_ind * sizeof(double), 
			   cudaMemcpyHostToDevice));

  cudaEventCreate(&start);
  cudaEventCreate(&stopKernel);
  cudaEventCreate(&stopMax);
  return 0;
}

//! @todo convert for gpu
void copyUpdateToDevice(unsigned id, unsigned iteration,  
			unsigned geno_count, unsigned n,
			char *d_snpMask, 
			int maxFIndex, double *d_Xtsnp, 
			size_t d_XtsnpPitch,
			const FortranMatrix &XtXi,
			const vector<double> &Xty){
  
  if(maxFIndex >= 0){
#ifdef _DEBUG
      cout << "iteration " << iteration << " id " << id 
	   << " masking index " << maxFIndex << " on device" 
	   << endl;
#endif
      cutilSafeCall(cudaMemsetAsync(d_snpMask + maxFIndex, 1, 1));
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
  /*! @todo fix for GPU
  cutilSafeCall(cudaMemcpy2D(d_Xtsnp + n, 
			     d_XtsnpPitch,
			     &XtSNP(n, (unsigned)0), 
			     (n + 1) * sizeof(double),
			     sizeof(double),
			     geno_count,
			     cudaMemcpyHostToDevice));
  */

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
