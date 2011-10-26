#include "type.h"

__device__ void reduceCore(const unsigned TID, unsigned N, double *reduce){
  /*
  if(!TID){
    for(int i = 1; i < N; i++)
      reduce[0] += reduce[i];
  }
  */
  ///*
  unsigned threads;
  while(N/2){
    if(N % 2){
      threads = (N + 1) / 2;
      //if(TID == threads - 1)
      //reduce[TID + threads] = 0;
      if(TID < threads - 1)
	reduce[TID] += reduce[TID + threads];
    } else {
      threads = N / 2;
      if(TID < threads)
	reduce[TID] += reduce[TID + threads];
    }
    __syncthreads();
    N = threads;
  }
  //*/
}

/*
  Written for register-based storage.  Result is placed in *reduce.
 */
__device__ void dotRR(const unsigned TID,
		    const unsigned N, 
		    const double x,
		    const double y,
		    double *reduce){

  // compute dot product terms
  // place x*y in shared memory (reduce)
  reduce[TID] = x * y;
  __syncthreads();
  reduceCore(TID, N, reduce);
}

__device__ void dotRG(const unsigned TID,
		      const unsigned N, 
		      const double x,
		      const double *y,
		      double *reduce){ // assume blockDim.x elements
  reduce[TID] = x * y[TID];
  __syncthreads();
  reduceCore(TID, N, reduce);
}

__device__ void dotGG(const unsigned TID,
		      const unsigned N, 
		      const double *x,
		      const double *y,
		      double *reduce){ // assume blockDim.x elements
  reduce[TID] = x[TID] * y[TID];
  __syncthreads();
  reduceCore(TID, N, reduce);
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
  for(int i = 0; i < N; i++){
    dotGG(TID, N, x, A + lda * i, reduce);
  }
  return reduce[TID];
}


__device__ double vecGMatG(const unsigned TID, 
			  const double *x,
			  const unsigned M, 
			  const unsigned N,
			  const double *A, 
			  const unsigned lda,
			  double *reduce) // reduce must be of length >= N, N <= M
{
  for(int i = 0; i < N; i++){
    dotRG(TID, M, x[TID], A + lda * i, reduce);
  }
  return reduce[TID];
}

