#include "type.h"

__device__ void reduceCore(const unsigned TID, unsigned N, ftype *reduce){
  unsigned threads;
  while(N){
    if(N % 2){
      threads = (N + 1) / 2;
      if(TID == threads - 1)
	reduce[TID + threads] = 0;
    } else {
      threads = N / 2;
    }
    if(TID < threads)
      reduce[TID] += reduce[TID + threads];
    __syncthreads();
    N /= 2;
  }
}

/*
  Written for register-based storage.  Result is placed in *reduce.
 */
__device__ void dotRR(const unsigned TID,
		    const unsigned N, 
		    const ftype x,
		    const ftype y,
		    ftype *reduce){

  // compute dot product terms
  // place x*y in shared memory (reduce)
  reduce[TID] = x * y;
  __syncthreads();
  reduceCore(TID, N, reduce);
}

__device__ void dotRG(const unsigned TID,
		      const unsigned N, 
		      const ftype x,
		      const ftype *y,
		      ftype *reduce){
  reduce[TID] = x * y[TID];
  __syncthreads();
  reduceCore(TID, N, reduce);
}


/*
  A is in constant memory.  A must be column-major.
 */
__device__ ftype vecRMatC(const unsigned TID,
			  const ftype x,
			  const unsigned M, 
			  const unsigned N,
			  const ftype *A, 
			  const unsigned lda,
			  ftype *reduce){
  ftype retVal;
  for(int i = 0; i < N; i++){
    dotRG(TID, M, x, A + lda * i, reduce);
    if(i == TID)
      retVal = *reduce;
  }
  return retVal;
}


__device__ ftype vecGMatG(const unsigned TID, 
			  const ftype *x,
			  const unsigned M, 
			  const unsigned N,
			  const ftype *A, 
			  const unsigned lda,
			  ftype *reduce){
  ftype retVal;
  for(int i = 0; i < N; i++){
    dotRG(TID, M, x[TID], A + lda * i, reduce);
    if(i == TID)
      retVal = *reduce;
  }
  return retVal;
}

