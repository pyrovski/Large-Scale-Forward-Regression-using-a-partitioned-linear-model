#include <cuda.h>
#include <cutil_inline.h>
#include <cublas.h>
#include <sstream>
#include <iostream>

#include "type.h"

#include "cuda_blas.cu"

#include "plm.h"

#define printBIDs 1

using namespace std;
//__shared__ ftype fval; // scalar
extern __shared__ ftype shared[];

__constant__ ftype d_Xty[iterationLimit + 27];

//const unsigned GPitch = (iterationLimit + 27) % 16 ? iterationLimit + 27 : (iterationLimit + 27)/16 * 16 + 16;

// extra space for padding.  Could make this triangular.
__constant__ ftype d_G[(iterationLimit + 27)*(iterationLimit + 27)];

__global__ void plm(// inputs
		    const unsigned m,          // rows of X
		    //const unsigned n,        // colums of X == number of blocks
		    const ftype *snptsnp,      // scalar, unique to block
		    const ftype *Xtsnp,        // n x 1 vector, unique to block
		    const unsigned XtsnpPitch, 
		    const ftype errorSS,       // scalar
		    const unsigned errorDF,       // scalar
		    //const ftype *G,          // symmetric matrix in const mem
		    //const ftype *Xty,        // n x 1 vector in const mem
		    const ftype *snpty,        // scalar, unique to block
		    const unsigned *snpMask,   // n x 1 vector
		    // outputs
		    ftype *f){
  /*! @todo could compute two SNPs per thread block.  
    This would ease the limitation of 8 thread blocks/MP for SM 1.3 devices,
    but might need some thread padding for warps.
   */

  ftype *reduce = shared; // n x 1
  //ftype *reduce2 = reduce + n;
  ftype GtXtsnp; // each thread stores one element of each array // Xtsnp
  //! @todo these might use fewer registers if kept in shared memory
  ftype snptmy; // scalar
  ftype s; // scalar

  unsigned BID = blockIdx.x + gridDim.x * blockIdx.y;
  unsigned TID = threadIdx.x;
  if(snpMask[BID]){
    // don't compute new F
    if(!TID)
      f[BID] = -1;
    return; 
  }
  // snptsnp - snptXGXtsnp:


  // GtXtsnp
  GtXtsnp = vecGMatCSq(TID, Xtsnp + BID * XtsnpPitch/sizeof(ftype), blockDim.x, d_G, 
		     blockDim.x,  //! length of column plus padding (no padding)
		     reduce); 
  
  // snptsnp - snptXGXtsnp
  dotRG(TID, blockDim.x, GtXtsnp, Xtsnp, reduce);
  s = snptsnp[BID] - *reduce;
#ifdef _DEBUG
  #if __CUDA_ARCH__ >= 200
  if(BID < printBIDs){
    printf("b%u\tt%u\tXtsnp: %le\n", BID, TID, Xtsnp[BID * XtsnpPitch/sizeof(ftype) + TID]);
    printf("b%u\tt%u\tGtXtsnp: %le\n", BID, TID, GtXtsnp);
    if(!TID){
      printf("b%u\tt%u\ts: %le\n", BID, TID, s);
    }
  }
  #endif
#endif
  // 1/(above)
  if(s > ftypeTol){
    s = (ftype)1/s;
    
    // snptmy
    dotRG(TID, blockDim.x, GtXtsnp, d_Xty, reduce);
    snptmy = -*reduce;
    
    if(!TID){
#ifdef _DEBUG
#if __CUDA_ARCH__ >= 200
      if(BID < printBIDs){
	printf("b%u\tt%u\tsnptXGXty: %le\n", BID, TID, snptmy);
	printf("b%u\tt%u\tsnpty: %le\n", BID, TID, snpty[BID]);
      }
#endif
#endif
      snptmy += snpty[BID];
      ftype modelSS = snptmy * snptmy * s;
      ftype errorSS2 = errorSS - modelSS;
      unsigned V2 = errorDF - 1;
      f[BID] = modelSS / errorSS2 * V2;
#ifdef _DEBUG
  #if __CUDA_ARCH__ >= 200
  if(BID < printBIDs){

    printf("b%u\tt%u\tmodelSS: %le\n", BID, TID, modelSS);
    printf("b%u\tt%u\tnew errorSS: %le\n", BID, TID, errorSS2);
    printf("b%u\tt%u\tnew V2: %u\n", BID, TID, V2);
    printf("b%u\tt%u\tf: %le\n", BID, TID, f[BID]);
  }
  #endif
#endif

    }
    return;
  } else {
    if(!TID)
      f[BID] = -1;
    return;
  }
}

cudaEvent_t start, stopKernel, stopMax;

unsigned plm_GPU(unsigned geno_count, unsigned blockSize, unsigned sharedSize, 
		 unsigned m, ftype* d_snptsnp, ftype* d_Xtsnp, 
		 unsigned d_XtsnpPitch, ftype ErrorSS, unsigned V2, 
		 ftype* d_snpty, unsigned* d_snpMask, ftype* d_f) throw(int)
{
    cublasGetError();
    cudaEventRecord(start, 0);
    plm<<<geno_count, blockSize, blockSize * sizeof(ftype)>>>
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
		 ftype *&d_snptsnp, ftype *&d_Xtsnp, size_t &d_XtsnpPitch, 
		 ftype *&d_snpty, unsigned *&d_snpMask, ftype *&d_f,
		 const vector<double> &SNPtSNP, const FortranMatrix &XtSNP,
		 const vector<double> &SNPty,
		 const vector<double> &Xty, const FortranMatrix &XtXi, 
		 const vector<unsigned> &snpMask){
  cudaEventCreate(&start);
  cudaEventCreate(&stopKernel);
  cudaEventCreate(&stopMax);

  cutilSafeCall(cudaMalloc(&d_snpMask, geno_count * sizeof(unsigned)));
  cutilSafeCall(cudaMemcpy(d_snpMask, &snpMask[0], 
			   geno_count * sizeof(unsigned), 
			   cudaMemcpyHostToDevice));

  
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

}

void copyUpdateToDevice(unsigned geno_count, unsigned n,
		       unsigned *d_snpMask, 
		       unsigned maxFIndex, ftype *d_Xtsnp, 
		       size_t d_XtsnpPitch,
		       const vector<unsigned> &snpMask,
		       FortranMatrix &XtSNP, const FortranMatrix &XtXi,
		       const vector<double> &Xty){
    cutilSafeCall(cudaMemcpy(d_snpMask + maxFIndex, &snpMask[maxFIndex], 
			     sizeof(unsigned), cudaMemcpyHostToDevice));

    //! copy updated XtSNP to GPU
    cutilSafeCall(cudaMemcpy2D(d_Xtsnp + n, 
			       d_XtsnpPitch,
			       &XtSNP(n, (unsigned)0), 
			       (n + 1) * sizeof(ftype),
			       sizeof(ftype),
			       geno_count,
			       cudaMemcpyHostToDevice));

    // update GPU G (const mem)
    cutilSafeCall(cudaMemcpyToSymbol(d_G, &XtXi.values[0], 
				     (n + 1) * (n + 1) * sizeof(ftype)));

    // update GPU Xty (const mem)
    // call fails unless we update the whole thing
    cutilSafeCall(cudaMemcpyToSymbol(d_Xty, &Xty[0], (n + 1) * sizeof(ftype)));

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

void getMaxF(unsigned id, unsigned iteration, unsigned geno_count, 
	     vector<double> &Fval, 
	     unsigned maxFIndex, ftype *d_f){
#ifndef _DEBUG
    cutilSafeCall(cudaMemcpy(&Fval[maxFIndex], &d_f[maxFIndex], sizeof(ftype),
			     cudaMemcpyDeviceToHost));
#else
    cutilSafeCall(cudaMemcpy(&Fval[0], d_f, geno_count * sizeof(ftype),
			     cudaMemcpyDeviceToHost));
    {
      stringstream ss;
      ss << "Fval_" << iteration << "_" << id << ".dat";
      writeD(ss.str(), Fval);
    }

#endif
}
