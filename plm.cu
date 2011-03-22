#include <cuda.h>
#include "type.h"

#include "cuda_blas.cu"

//__shared__ ftype fval; // scalar
extern __shared__ ftype shared[];

uint shared_size(uint n){
  return n;
}
__constant__ ftype d_Xty[iterationLimit + 27];

//const unsigned GPitch = (iterationLimit + 27) % 16 ? iterationLimit + 27 : (iterationLimit + 27)/16 * 16 + 16;

// extra space for padding.  Could make this triangular.
__constant__ ftype d_G[(iterationLimit + 27)*(iterationLimit + 27)];

//! @todo could replace N with blockDim.x
__global__ void plm(// inputs
		    const unsigned m,          // rows of X
		    //const unsigned n,        // colums of X
		    //const ftype *X,            // m x n matrix. padding?
		    //const ftype *snp,        // m x 1 vector, unique to block
		    //const unsigned XPitch,   // 
		    const ftype *snptsnp,      // scalar, unique to block
		    const ftype *Xtsnp,        // n x 1 vector, unique to block
		    const unsigned XtsnpPitch, 
		    const ftype errorSS,       // scalar
		    const ftype errorDF,       // scalar
		    //const ftype *G,          // symmetric matrix. padding?
		    //const ftype *Xty,        // n x 1 vector
		    const ftype *snpty,        // scalar, unique to block
		    // outputs
		    ftype *f){
  ftype *reduce = shared; // n x 1
  //ftype *reduce2 = reduce + n;
  ftype GtXtsnp; // each thread stores one element of each array // Xtsnp
  //! @todo these might use fewer registers if kept in shared memory
  ftype snptmy; // scalar
  ftype s; // scalar

  unsigned BID = blockIdx.x + gridDim.x * blockIdx.y;
  unsigned TID = threadIdx.x;
  // snptsnp - snptXGXtsnp:

  // (snptX)' = Xtsnp

  // Xtsnp
  /*! @todo maintain this on the host
     @todo as written, this needs m elements in shared memory, not n
  Xtsnp = vecGMatG(TID, 
		   snp + BID * snpPitch,
		   m, n, X, 
		   snpPitch,   //! length of column plus padding
		   reduce); 
  
  */
  // GtXtsnp
  GtXtsnp = vecGMatCSq(TID, Xtsnp, blockDim.x, d_G, 
		     blockDim.x,  //! length of column plus padding
		     reduce); 

  // snptsnp - snptXGXtsnp
  dotRG(TID, blockDim.x, GtXtsnp, Xtsnp, reduce);
  s = snptsnp[BID] - *reduce;
  // 1/(above)
  if(s > ftypeTol){
    s = (ftype)1/s;
    
    // snptmy
    dotRG(TID, blockDim.x, GtXtsnp, d_Xty, reduce);
    snptmy = *reduce;
    
    if(!TID){
      snptmy += snpty[BID];
      ftype modelSS = snptmy * snptmy * s;
      ftype errorSS2 = errorSS - modelSS;
      ftype V2 = errorDF - 1;
      f[BID] = modelSS / errorSS2 * V2;
    }
    return;
  } else {
    if(!TID)
      f[BID] = 0;
    return;
  }
}
