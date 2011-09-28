#include <vector>

#include "type.h"
#include "fortran_matrix.h"

/*! @todo make hardIterationLimit a hard limit, and iterationLimit a soft limit.
  hardIterationLimit is constrained by constant memory size, which varies
  according to device family, so we can define hardIterationLimit as 
  a macro function of 
 */
#define iterationLimit 50

unsigned plm_GPU(unsigned geno_count, unsigned blockSize, 
		 unsigned m, ftype* d_snptsnp, ftype* d_Xtsnp, 
		 unsigned d_XtsnpPitch, ftype ErrorSS, unsigned V2, 
		 ftype* d_snpty, unsigned* d_snpMask, ftype* d_f,
		 std::vector<double> &Fval) throw(int);

int copyToDevice(const unsigned id, 
		 const unsigned verbosity,
		 const unsigned geno_count, const unsigned n,
		 ftype *&d_snptsnp, ftype *&d_Xtsnp, size_t &d_XtsnpPitch, 
		 ftype *&d_snpty, unsigned *&d_snpMask, ftype *&d_f,
		  const std::vector<double> &SNPtSNP, const FortranMatrix &XtSNP,
		 const std::vector<double> &SNPty,
		 const std::vector<double> &Xty, const FortranMatrix &XtXi, 
		 const std::vector<unsigned> &snpMask);

void copyUpdateToDevice(unsigned id, unsigned iteration,
			unsigned geno_count, unsigned n,
			unsigned *d_snpMask, 
			int maxFIndex, ftype *d_Xtsnp, 
			size_t d_XtsnpPitch,
			const std::vector<unsigned> &snpMask,
			FortranMatrix &XtSNP, const FortranMatrix &XtXi,
			const std::vector<double> &Xty);

float getGPUCompTime();

float getGPUMaxTime();

void getMaxFGPU(unsigned id, unsigned iteration, unsigned geno_count, 
	     std::vector<double> &Fval, 
	     unsigned maxFIndex, ftype *d_f);
