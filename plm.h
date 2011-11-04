#include <vector>

#include "type.h"
#include "fortran_matrix.h"

/*! @todo make hardIterationLimit a hard limit, and iterationLimit a soft limit.
  hardIterationLimit is constrained by constant memory size, which varies
  according to device family, so we can define hardIterationLimit as 
  a macro function of 
 */
#define fixedPlusIteration_limit 89

unsigned plm_GPU(unsigned geno_count, unsigned blockSize, 
		 unsigned m, double* d_snptsnp, double* d_Xtsnp, 
		 unsigned d_XtsnpPitch, double ErrorSS, unsigned V2, 
		 double* d_snpty, char* d_snpMask, float* d_f,
		 std::vector<float> &Fval) throw(int);

int copyToDevice(const unsigned id, 
		 const unsigned verbosity,
		 const unsigned geno_count, const unsigned n,
		 double *&d_snptsnp, double *&d_Xtsnp, size_t &d_XtsnpPitch, 
		 double *&d_snpty, char *&d_snpMask, float *&d_f,
		  const std::vector<double> &SNPtSNP, const FortranMatrix &XtSNP,
		 const std::vector<double> &SNPty,
		 const std::vector<double> &Xty, const FortranMatrix &XtXi, 
		 const std::vector<char> &snpMask);

void copyUpdateToDevice(unsigned id, unsigned iteration,
			unsigned geno_count, unsigned n,
			char *d_snpMask, 
			int maxFIndex, double *d_Xtsnp, 
			size_t d_XtsnpPitch,
			const std::vector<char> &snpMask,
			FortranMatrix &XtSNP, const FortranMatrix &XtXi,
			const std::vector<double> &Xty);

float getGPUCompTime();

float getGPUMaxTime();

void getMaxFGPU(unsigned id, unsigned iteration, unsigned geno_count, 
	     std::vector<float> &Fval, 
	     unsigned maxFIndex, float *d_f);
