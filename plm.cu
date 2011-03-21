#include <cuda.h>
#include "type.h"

#include "cuda_blas.cu"

extern __shared__ ftype shared[];

__global__ void plm(// inputs
		    const unsigned n,    // colums of X
		    const unsigned m,    // rows of X
		    const ftype *X,      // m x n matrix. padding?
		    const ftype *snp,    // m x 1 vector, unique to block
		    const ftype *snptsnp, // scalar, unique to block
		    const ftype errorSS, // scalar
		    const ftype errorDF, // scalar
		    const ftype *G,      // symmetric matrix. padding?
		    const ftype *Xty,    // n x 1 vector
		    const ftype *snpty,   // scalar, unique to block
		    // outputs
		    ftype *f){

  ftype *Xtsnp = shared; // n x 1
  ftype *GtXtsnp = Xtsnp + n; // n x 1.  G is symmetric, so doesn't need a transpose
  ftype *s = GtXtsnp + n; // scalar
  ftype *fval = s + 1; // scalar
  ftype *snptmy = fval + 1; // scalar

  unsigned BID = blockIdx.x + gridDim.x * blockIdx.y;
  unsigned TID = threadIdx.x;
  // snptsnp - snptXGXtsnp:

  // (snptX)' = Xtsnp
  // don't worry about aligning first thread to 128-byte boundary; assume compute capability 1.2+
  matVec(m, n, X, 
	 m,  //! @todo length of column plus padding
	 snp + BID * m, 
	 Xtsnp); 

  // GtXtsnp
  matVec(n, n, G, 
	 n,  //! @todo length of column plus padding
	 Xtsnp, GtXtsnp); 

  // snptsnp - snptXGXtsnp
  dot(n, Xtsnp, GtXtsnp, s);
  if(!TID)
    *s = snptsnp[BID] - *s;
  __syncthreads();
  // 1/(above)
  if(*s > ftypeTol){
    if(!TID)
      *s = 1/ *s;
    __syncthreads();

    dot(n, GtXtsnp, Xty, snptmy);
    
    if(!TID){
      *snptmy += snpty[BID];
      ftype modelSS = *snptmy * *snptmy * *s;
      ftype errorSS2 = errorSS - modelSS;
      ftype V2 = errorDF - 1;
      *fval = modelSS / errorSS2 * V2;
    }
    __syncthreads();
  } else {
    if(!TID)
      fval = 0;
    __syncthreads();
  }
  f[BID] = *fval;
}
