#include <cuda.h>
#include "type.h"

#include "cuda_blas.cu"

//__shared__ ftype fval; // scalar
extern __shared__ ftype shared[];

__constant__ ftype d_Xty[iterationLimit + 27];

//const unsigned GPitch = (iterationLimit + 27) % 16 ? iterationLimit + 27 : (iterationLimit + 27)/16 * 16 + 16;

// extra space for padding.  Could make this triangular.
__constant__ ftype d_G[(iterationLimit + 27)*(iterationLimit + 27)];

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
		    const unsigned errorDF,       // scalar
		    //const ftype *G,          // symmetric matrix. padding?
		    //const ftype *Xty,        // n x 1 vector
		    const ftype *snpty,        // scalar, unique to block
		    const unsigned *snpMask,
		    // outputs
		    ftype *f){
  /*! @todo could compute two SNPs per thread block.  
    This would ease the limitation of 8 thread blocks/MP for SM 1.3 devices.
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
  GtXtsnp = vecGMatCSq(TID, Xtsnp, blockDim.x, d_G, 
		     blockDim.x,  //! length of column plus padding (no padding)
		     reduce); 
  
  // snptsnp - snptXGXtsnp
  dotRG(TID, blockDim.x, GtXtsnp, Xtsnp, reduce);
  s = snptsnp[BID] - *reduce;
#ifdef _DEBUG
  #if __CUDA_ARCH__ >= 200
  if(!BID){
    printf("b%u\tt%u\tGtXtsnp: %le\n", BID, TID, GtXtsnp);
    printf("b%u\tt%u\tXtsnp: %le\n", BID, TID, Xtsnp[TID]);
    printf("b%u\tt%u\tG(1,%u): %le\n", BID, TID, TID, d_G[TID]);
    if(!TID){
      printf("b%u\tt%u\tsnptsnp: %le\n", BID, TID, snptsnp[BID]);
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
      f[BID] = -1;
    return;
  }
}
