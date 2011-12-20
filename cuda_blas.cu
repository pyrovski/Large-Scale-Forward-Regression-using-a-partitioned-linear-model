#include "type.h"

extern __shared__ double shared[];
#if __CUDA_ARCH__ >= 200
#define sharedSize 49152
#else
#define sharedSize 16384
#endif

template <unsigned blockSize> __device__ 
void warpReduce(const unsigned TID, volatile double *reduce){
  // warp-synchronous, with max 32 threads!
  if(blockSize >= 64) reduce[TID] += reduce[TID + 32];
  if(blockSize >= 32) reduce[TID] += reduce[TID + 16];
  if(blockSize >= 16) reduce[TID] += reduce[TID + 8];
  if(blockSize >= 8)  reduce[TID] += reduce[TID + 4];
  if(blockSize >= 4)  reduce[TID] += reduce[TID + 2];
  if(blockSize >= 2)  reduce[TID] += reduce[TID + 1];
}

template <unsigned blockSize> 
__device__ void reduceCorePow2(const unsigned TID, double *reduce){
  if(blockSize == 512){if(TID < 256){reduce[TID] += reduce[TID + 256];} __syncthreads();}
  if(blockSize >= 256){if(TID < 128){reduce[TID] += reduce[TID + 128];} __syncthreads();}
  if(blockSize >= 128){if(TID < 64) {reduce[TID] += reduce[TID +  64];} __syncthreads();}
  if(TID < 32) warpReduce<blockSize>(TID, reduce);
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
  unsigned threads;
  while(N/2){
    threads = (N + 1) / 2;
    //! @todo can remove this condition test if threads <= 64, 
    //as it doesn't affect the result or save time,
    //assuming shared size is 2*threads
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
		    const double y){

  double *reduce = shared;

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
		      const double *y){
  double *reduce = shared;
  __syncthreads();
  reduce[TID] = x * y[TID];
  __syncthreads();
  reduceCore(TID, N, reduce);
}

/*
  A is in constant memory.  A must be column-major and square (NxN).
*/
__device__ double vecRMatCSq(const unsigned TID,
			     const double x,
			     const unsigned N,
			     const double *A, 
			     const unsigned lda){
    
  double retVal = 0.0;
  __syncthreads();
  shared[TID] = x;
  __syncthreads();
  /*
    for(int i = 0; i < N; i++){
    if(i == TID)
      continue;
    retVal += shared[i] * A[lda * TID + i];
  }
  retVal += x * A[lda * TID + TID];
  */
  for(int i = 0; i < N; i++)
    retVal += shared[i] * A[lda * i + TID];
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
template<unsigned blockSize> __global__ void 
columnDot(const double *d_mat, unsigned n_rows, unsigned n_cols, 
	  unsigned columnPitchInWords, // words between tops of columns
	  double *result, 
	  unsigned resultStrideInWords // words between result elements
	  ){
  __shared__ double myShared[blockSize];
  const unsigned BID = blockIdx.x + gridDim.x * blockIdx.y;
  const unsigned TID = threadIdx.x;
  
  const unsigned col = BID;
  double myResult = 0.0;

  for(unsigned row = 0; row < n_rows; row += blockSize)
    if(row + TID < n_rows)
      myResult += d_mat[col * columnPitchInWords + row + TID] * d_mat[col * columnPitchInWords + row + TID];

  myShared[TID] = myResult;
  __syncthreads();
  reduceCorePow2<blockSize>(TID, myShared);
  if(!TID)
    result[BID * resultStrideInWords] = myShared[0];
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
    std::cerr << "invalid blockSize:" << blockSize << std::endl;
    exit(1);
  }

  cudaEventRecord(start, 0);
  columnDot<blockSize>
    <<<grid, blockSize>>>
    (d_mat, n_rows, n_cols, columnPitchInWords, d_result, resultStrideInWords);
  cudaEventRecord(stopKernel, 0);

  #warning synchronizing for timing
  cudaEventSynchronize(stopKernel);
}
