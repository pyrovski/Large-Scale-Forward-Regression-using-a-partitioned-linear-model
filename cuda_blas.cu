#include "type.h"

extern __shared__ double shared[];
#if __CUDA_ARCH__ >= 200
#define sharedSize 49152
#else
#define sharedSize 16384
#endif

template <unsigned blockSize> __device__ 
void warpReduce(unsigned TID, volatile double *reduce){
  // warp-synchronous, with max 32 threads!
  if(blockSize >= 64) reduce[TID] += reduce[TID + 32];
  if(blockSize >= 32) reduce[TID] += reduce[TID + 16];
  if(blockSize >= 16) reduce[TID] += reduce[TID + 8];
  if(blockSize >= 8)  reduce[TID] += reduce[TID + 4];
  if(blockSize >= 4)  reduce[TID] += reduce[TID + 2];
  if(blockSize >= 2)  reduce[TID] += reduce[TID + 1];
}

template <unsigned blockSize> __device__ void reduceCorePow2(const unsigned TID, double *reduce){
  if(blockSize == 512){if(TID < 256){reduce[TID] += reduce[TID + 256];} __syncthreads();}
  if(blockSize >= 256){if(TID < 128){reduce[TID] += reduce[TID + 128];} __syncthreads();}
  if(blockSize >= 127){if(TID < 64) {reduce[TID] += reduce[TID +  64];} __syncthreads();}
  if(TID < 32) warpReduce(TID, reduce);
}

/*
__device__ void reduceMadDotGG(const unsigned TID, unsigned N, double *reduce, 
			  const double *x, const double *y){
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
*/

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
    /*
    if(N % 2){
      threads = (N + 1) / 2;
      //if(TID == threads - 1)
      //reduce[TID + threads] = 0;
      //! @todo can remove this condition test if threads <= 64, 
      //as it doesn't affect the result or save time,
      //assuming shared size is 2*threads

      if(TID < threads - 1) 
	reduce[TID] += reduce[TID + threads];
    } else {
      threads = N / 2;
      if(TID < threads)
	reduce[TID] += reduce[TID + threads];
    }
    */
    threads = (N + 1) / 2;
    if(TID < threads - (N % 2))
      reduce[TID] += reduce[TID + threads];
    
    __syncthreads();
    N = threads;
  }
  //*/
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
__device__ double vecRMatCSq(const unsigned TID,
			  const double x,
			  const unsigned N,
			  const double *A, 
			  const unsigned lda,
			  double *reduce){
  double retVal;
  for(int i = 0; i < N; i++){
    dotRG(TID, N, x, A + lda * i, reduce);
    if(i == TID)
      retVal = *reduce;
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


/* ignore this; it is a mess...
__global__ template<unsigned blockSize> void 
columnDot___(const double *d_mat, unsigned n_rows, uint64_t n_cols, 
	  unsigned columnPitchInWords, // words between tops of columns
	  double *result, 
	  unsigned resultStrideInWords // words between result elements
	  ){
  const unsigned BID = blockIdx.x + gridDim.x * blockIdx.y;
  const unsigned TID = threadIdx.x;
  
  // if GID >= n_cols, participate in reading, but not calculating or writing
  const unsigned GID = BID * blockSize + TID;
  const unsigned m = sharedSize / sizeof(double) / blockSize;

  unsigned col = BID;
  double myResult = 0.0;

#ifdef splitDoubles
#else
  for(unsigned row = 0; row < n_rows; row += m){
    for(unsigned colOffset = 0; col + colOffset < min(col + blockSize, n_cols); colOffset++){
      if(TID < m && row + TID < n_rows){
#ifdef shared1
	// bank conflicts on write
	shared[TID * m + colOffset] = 
          d_mat[(col + colOffset)*columnPitchInWords + row + TID];
#else
	// minimal conflicts on write
	shared[TID + m * colOffset] = 
	  d_mat[(col + colOffset)*columnPitchInWords + row + TID];
#endif
      }
    }
    __syncthreads();
#ifdef shared1
    for(unsigned rowOffset = 0; row + rowOffset < min(n_rows, row + m); rowOffset++){
      // minimal bank conflicts on read
      myResult += shared[TID + m * rowOffset] * shared[TID + m * rowOffset];
    }
#else
    // bank conflicts on read
    //for(unsigned rowOffset = 0; row + rowOffset < min(n_rows, row + m); rowOffset++)
    //myResult += shared[TID * m + rowOffset] * shared[TID * m + rowOffset];
    
    // minimal bank conflicts on read
    unsigned maxRowCount = min(n_rows, row + m) - row;
    for(unsigned rowOffset = TID % m, rowCount = 0; rowCount < maxRowCount ; rowCount++){
      myResult += shared[rowOffset + TID * m] + shared[rowOffset + TID * m];
      rowOffset = (rowOffset + 1) % m;
    }
#endif
    __syncthreads();
  }
#endif // splitDoubles
  result[GID * resultStrideInWords] = myResult;
}
*/

/*
  Compute the dot product of each column with itself.
  Matrix is stored in column-major order.
 */
__global__ template<unsigned blockSize> void 
columnDot(const double *d_mat, unsigned n_rows, unsigned n_cols, 
	  unsigned columnPitchInWords, // words between tops of columns
	  double *result, 
	  unsigned resultStrideInWords // words between result elements
	  ){
  const unsigned BID = blockIdx.x + gridDim.x * blockIdx.y;
  const unsigned TID = threadIdx.x;
  
  const unsigned col = BID;
  double myResult = 0.0, tmp;

  for(unsigned row = 0; row < n_rows; row += blockSize){
    if(row + TID < n_rows){
      tmp = 0.0;
    }else{
      tmp = d_mat[col * columnPitchInWords + row + TID];
    }
    myResult += tmp * tmp;
  }

  shared[TID] = myResult;
  __syncthreads();
  reduceCorePow2<blockSize>(TID, blockSize, shared);
  if(!TID)
    result[BID * resultStrideInWords] = shared[0];
}

void columnDot_gpu(const double *d_mat, unsigned n_rows, uint64_t n_cols, 
		   unsigned columnPitchInWords, // words between tops of columns
		   double *d_result, 
		   unsigned resultStrideInWords // words between result elements
		   ){
  dim3 grid;
  initGrid(grid, n_cols);
  const unsigned blockSize = 256;

  if(blockSize <= 32 || blockSize > 512){
    cerr << "invalid blockSize:" << blockSize << endl;
    exit(1);
  }
  columnDot<blockSize>
    <<<grid, blockSize, blockSize * sizeof(double)>>>
    (d_mat, n_rows, n_cols, columnPitchInWords, d_result, resultStrideInWords);
}
