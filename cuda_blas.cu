#include "type.h"

__device__ void reduceCorePow2(const unsigned TID, unsigned N, double *reduce){
}

__device__ void reduceMadDotGG(const unsigned TID, unsigned N, double *reduce, 
			  const double *x, const double *y){
  // serial version for clarity
  /*
  if(!TID){
    reduce[0] = x[0] * y[0];
      for(int i = 1; i < N; i++)
	reduce[0] += x[i]*y[i];
  }
  */

  unsigned threads, nextSmaller;
  // parallel version
  if(__popc(N) != 1){ // non power of two
    // do one reduction to a power of two
    nextSmaller = 1 << (31-__clz(N));
    threads = nextSmaller/2;
    //! @todo threads could be > nextSmaller/2,
    if(TID < threads){
      reduce[TID] = x[TID] * y[TID] + x[TID + threads] * y[TID + threads];
      if(TID < N - nextSmaller)
	reduce[TID] += x[TID + nextSmaller] * y[TID + nextSmaller];
    }
    __syncthreads();
  } else {
    threads = N/2;
  }
  reduceCorePow2(TID, threads, reduce);
}


__device__ void reduceCore(const unsigned TID, unsigned N, double *reduce){
  /*
  if(!TID){
    for(int i = 1; i < N; i++)
      reduce[0] += reduce[i];
  }
  */
  unsigned threads;
  while(N/2){
    threads = (N + 1) / 2;
    if(TID < threads - (N % 2))
      reduce[TID] += reduce[TID + threads];
    
    __syncthreads();
    N = threads;
  }
}

/*
  Written for register-based storage.  Result is placed in *reduce.
 */
/*! @todo dot products could be faster; specifically, each thread is only 
  multiplying or adding, not both simultaneously.
 */
__device__ void dotRR(const unsigned TID,
		    const unsigned N, 
		    const double x,
		    const double y,
		    double *reduce){

  // compute dot product terms
  // place x*y in shared memory (reduce)
  __syncthreads();
  reduce[TID] = x * y;
  __syncthreads();
  reduceCore(TID, N, reduce);
}

__device__ void dotRG(const unsigned TID,
		      const unsigned N, 
		      const double x,
		      const double *y,
		      double *reduce){ // assume blockDim.x elements
  __syncthreads();
  reduce[TID] = x * y[TID];
  __syncthreads();
  reduceCore(TID, N, reduce);
}

__device__ void dotGG(const unsigned TID,
		      const unsigned N, 
		      const double *x,
		      const double *y,
		      double *reduce){ // assume blockDim.x elements
  __syncthreads();
  reduce[TID] = x[TID] * y[TID];
  __syncthreads();
  reduceCore(TID, N, reduce);
}

/*
  A is in constant memory.  A must be column-major and square (NxN).
 */
__device__ double vecRMatCSq(const unsigned TID,
			     const unsigned BID,
			  const double x,
			  const unsigned N,
			  const double *A, 
			  const unsigned lda,
			  double *reduce){
  double retVal = 0.0;
  __syncthreads();
  reduce[TID] = x;
  __syncthreads();
  /*
  for(int i = 0; i < N; i++){
    if(i == TID)
      continue;
    retVal += reduce[i] * A[lda * TID + i];
  }
  retVal += x * A[lda * TID + TID];
  */
  for(int i = 0; i < N; i++){
    retVal += reduce[i] * A[lda * TID + i];
    #ifdef printGPU
    if(printBIDs(BID)){
      if(!TID)
	printf("b%03u\tt%03u\tmultiplying:\t%1.10le\t*\t%1.10le:\t%1.10le,\tsum:%1.10le\n", BID, TID, reduce[i], A[lda * TID + i], reduce[i] * A[lda * TID + i], retVal);
    }
    #endif
  }
  return retVal;
}

/*
  A is in constant memory.  A must be column-major and square.
 */
__device__ double vecGMatCSq(const unsigned TID,
			  const double *x,
			  const unsigned N,
			  const double *A, 
			  const unsigned lda,
			  double *reduce){
  double retVal;
  for(int i = 0; i < N; i++){
    dotGG(TID, N, x, A + lda * i, reduce);
    if(i == TID)
      retVal = *reduce;
  }
  return retVal;
}


__device__ double vecGMatG(const unsigned TID, 
			  const double *x,
			  const unsigned M, 
			  const unsigned N,
			  const double *A, 
			  const unsigned lda,
			  double *reduce) // reduce must be of length >= N, N <= M
{
  double retVal;
  for(int i = 0; i < N; i++){
    dotRG(TID, M, x[TID], A + lda * i, reduce);
    if(i == TID)
      retVal = *reduce;
  }
  return retVal;
}

