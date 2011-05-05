#include <cuda.h>
#include <cutil_inline.h>
#include <cublas.h>
#include <sstream>
#include <iostream>

#include "type.h"

#include "plm.h"

#ifdef _DEBUG
#if __CUDA_ARCH__ >= 200
//#define printGPU
#endif
#endif

__device__ int printBIDs(unsigned BID){
  return(BID == 0 || BID == 40);
}

#include "cuda_blas.cu"

using namespace std;
//__shared__ double fval; // scalar
extern __shared__ double shared[];

__constant__ double d_Xty[iterationLimit + 27];

//const unsigned GPitch = (iterationLimit + 27) % 16 ? iterationLimit + 27 : (iterationLimit + 27)/16 * 16 + 16;

// extra space for padding.  Could make this triangular.
__constant__ double d_G[(iterationLimit + 27)*(iterationLimit + 27)];

__global__ void plm(// inputs
		    const unsigned m,          // rows of X
		    //const unsigned n,        // colums of X == number of blocks
		    const double *snptsnp,      // scalar, unique to block
		    const double *Xtsnp,        // n x 1 vector, unique to block
		    const unsigned XtsnpPitch, 
		    const double errorSS,       // scalar
		    const unsigned errorDF,       // scalar
		    //const double *G,          // symmetric matrix in const mem
		    //const double *Xty,        // n x 1 vector in const mem
		    const double *snpty,        // scalar, unique to block
		    const unsigned *snpMask,   // n x 1 vector
		    // outputs
		    double *f){
  /*! @todo could compute two SNPs per thread block.  
    This would ease the limitation of 8 thread blocks/MP for SM 1.3 devices,
    but might need some thread padding for warps.
   */

  double *reduce = shared; // n x 1
  //double *reduce2 = reduce + n;
  double GtXtsnp; // each thread stores one element of each array // Xtsnp
  //! @todo these might use fewer registers if kept in shared memory
  double snptmy; // scalar
  double s; // scalar

  unsigned BID = blockIdx.x + gridDim.x * blockIdx.y;
  unsigned TID = threadIdx.x;
  if(snpMask[BID]){
    // don't compute new F
    if(!TID)
      f[BID] = 0;
    return; 
  }
  // snptsnp - snptXGXtsnp:


  // GtXtsnp
  GtXtsnp = vecGMatCSq(TID, Xtsnp + BID * XtsnpPitch/sizeof(double), blockDim.x, d_G, 
		     blockDim.x,  //! length of column plus padding (no padding)
		     reduce); 
  
  // snptsnp - snptXGXtsnp
  dotRG(TID, blockDim.x, GtXtsnp, Xtsnp + BID * XtsnpPitch/sizeof(double), reduce);
  s = snptsnp[BID] - *reduce;
#ifdef printGPU
  if(printBIDs(BID)){
    printf("b%03u\tt%03u\tXtsnp: %1.10le\n", BID, TID, Xtsnp[BID * XtsnpPitch/sizeof(double) + TID]);
    printf("b%03u\tt%03u\tGtXtsnp: %1.10le\n", BID, TID, GtXtsnp);
    if(!TID){
      printf("b%03u\tt%03u\tsnptsnp: %1.10le\n", BID, TID, snptsnp[BID]);
      printf("b%03u\tt%03u\tsnptXGXtsnp: %1.10le\n", BID, TID, *reduce);
      printf("b%03u\tt%03u\ts: %1.10le\n", BID, TID, s);
    }
  }
#endif
  // 1/(above)
  if(s > doubleTol){
    s = (double)1/s;
    
    // snptmy
    dotRG(TID, blockDim.x, GtXtsnp, d_Xty, reduce);
    snptmy = -*reduce;
    
    if(!TID){
#ifdef printGPU
      if(printBIDs(BID)){
	printf("b%03u\tt%03u\tsnptXGXty: %1.10le\n", BID, TID, snptmy);
	printf("b%03u\tt%03u\tsnpty: %1.10le\n", BID, TID, snpty[BID]);
      }
#endif
      snptmy += snpty[BID];
      double modelSS = snptmy * snptmy * s;
      double errorSS2 = errorSS - modelSS;
      unsigned V2 = errorDF - 1;
      f[BID] = modelSS / errorSS2 * V2;
#ifdef printGPU
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
}

cudaEvent_t start, stopKernel, stopMax, startCopy, stopCopy, 
  startCopyUpdate, stopCopyUpdate;

unsigned plm_GPU(unsigned geno_count, unsigned blockSize, 
		 unsigned m, double* d_snptsnp, double* d_Xtsnp, 
		 unsigned d_XtsnpPitch, double ErrorSS, unsigned V2, 
		 double* d_snpty, unsigned* d_snpMask, double* d_f) throw(int)
{
    cublasGetError();
    cudaEventRecord(start, 0);
    plm<<<geno_count, blockSize, blockSize * sizeof(double)>>>
      (m ,        
       d_snptsnp, 
       d_Xtsnp, 
       d_XtsnpPitch, 
       ErrorSS, V2, 
       d_snpty, 
       d_snpMask,
       d_f);
    cudaEventRecord(stopKernel, 0);
    cutilSafeCall(cudaThreadSynchronize());
    // cublas uses 1-based index
    unsigned maxFIndex = cublasIdamax(geno_count, d_f, 1);
    cudaEventRecord(stopMax, 0);
    if(!maxFIndex){
      cerr << "maxFIndex <= 0!" << endl;
      throw(1);
    }
    maxFIndex -= 1;

    return maxFIndex;
}

/*!
  should only be called once
 */
void copyToDevice(unsigned geno_count, const unsigned n, 
		 double *&d_snptsnp, double *&d_Xtsnp, size_t &d_XtsnpPitch, 
		 double *&d_snpty, unsigned *&d_snpMask, double *&d_f,
		 const vector<double> &SNPtSNP, const FortranMatrix &XtSNP,
		 const vector<double> &SNPty,
		 const vector<double> &Xty, const FortranMatrix &XtXi, 
		 const vector<unsigned> &snpMask){
  cudaEventCreate(&start);
  cudaEventCreate(&stopKernel);
  cudaEventCreate(&stopMax);
  cudaEventCreate(&startCopy);
  cudaEventCreate(&stopCopy);
  cudaEventCreate(&startCopyUpdate);
  cudaEventCreate(&stopCopyUpdate);

  cudaEventRecord(startCopy, 0);
  cutilSafeCall(cudaMalloc(&d_snpMask, geno_count * sizeof(unsigned)));
  cutilSafeCall(cudaMemcpy(d_snpMask, &snpMask[0], 
			   geno_count * sizeof(unsigned), 
			   cudaMemcpyHostToDevice));

  
  //! @todo this won't be coalesced
  cutilSafeCall(cudaMalloc(&d_snptsnp, geno_count * sizeof(double)));
  cutilSafeCall(cudaMemcpy(d_snptsnp, &SNPtSNP[0], geno_count * sizeof(double), 
			   cudaMemcpyHostToDevice));

  cutilSafeCall(cudaMallocPitch(&d_Xtsnp, &d_XtsnpPitch, 
				(n + iterationLimit) * sizeof(double), 
				geno_count));
  cutilSafeCall(cudaMemcpy2D(d_Xtsnp, d_XtsnpPitch, &XtSNP.values[0], 
			     n * sizeof(double), n * sizeof(double), geno_count, 
			     cudaMemcpyHostToDevice));
  
  cutilSafeCall(cudaMalloc(&d_snpty, geno_count * sizeof(double)));
  cutilSafeCall(cudaMemcpy(d_snpty, &SNPty[0], geno_count * sizeof(double), 
			   cudaMemcpyHostToDevice));
  
  cutilSafeCall(cudaMemcpyToSymbol(d_G, &XtXi.values[0], n * n * sizeof(double)));
  cutilSafeCall(cudaMemcpyToSymbol(d_Xty, &Xty[0], n * sizeof(double)));
  
  cutilSafeCall(cudaMalloc(&d_f, geno_count * sizeof(double)));
  cudaEventRecord(stopCopy, 0);
}

void copyUpdateToDevice(unsigned id, unsigned iteration,  
			unsigned geno_count, unsigned n,
		       unsigned *d_snpMask, 
		       int maxFIndex, double *d_Xtsnp, 
		       size_t d_XtsnpPitch,
		       const vector<unsigned> &snpMask,
		       FortranMatrix &XtSNP, const FortranMatrix &XtXi,
		       const vector<double> &Xty){

  cudaEventRecord(startCopyUpdate, 0);

  if(maxFIndex >= 0){
#ifdef _DEBUG
      cout << "iteration " << iteration << " id " << id 
	   << " masking index " << maxFIndex << " on device: " 
	   << snpMask[maxFIndex] << endl;
#endif
    cutilSafeCall(cudaMemcpy(d_snpMask + maxFIndex, &snpMask[maxFIndex], 
			     sizeof(unsigned), cudaMemcpyHostToDevice));
    unsigned maskVal;
    cutilSafeCall(cudaMemcpy(&maskVal, d_snpMask + maxFIndex,
			     sizeof(unsigned), cudaMemcpyDeviceToHost));
#ifdef _DEBUG
      cout << "iteration " << iteration << " id " << id 
	   << " mask index " << maxFIndex << ": "
	   << maskVal << endl;
#endif
      /* this doesn't help...
    cutilSafeCall(cudaMemcpy(d_snpMask, &snpMask[0], 
			     sizeof(unsigned) * geno_count, cudaMemcpyHostToDevice));
      */
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

    cudaEventRecord(stopCopyUpdate, 0);
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

float getGPUCopyTime(){
  float copy_elapsed_time;
  cudaEventElapsedTime(&copy_elapsed_time, startCopy, stopCopy);
  return copy_elapsed_time / 1000.0f;
}

float getGPUCopyUpdateTime(){
  float copy_elapsed_time;
  cudaEventElapsedTime(&copy_elapsed_time, startCopyUpdate, stopCopyUpdate);
  return copy_elapsed_time / 1000.0f;
}

void getMaxF(unsigned id, unsigned iteration, unsigned geno_count, 
	     vector<double> &Fval, 
	     unsigned maxFIndex, double *d_f){
#ifndef _DEBUG
    cutilSafeCall(cudaMemcpy(&Fval[maxFIndex], &d_f[maxFIndex], sizeof(double),
			     cudaMemcpyDeviceToHost));
#else
    cutilSafeCall(cudaMemcpy(&Fval[0], d_f, geno_count * sizeof(double),
			     cudaMemcpyDeviceToHost));
    {
      stringstream ss;
      ss << "Fval_" << iteration << "_" << id << ".dat";
      writeD(ss.str(), Fval);
    }

#endif
}
