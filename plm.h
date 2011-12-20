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

int copyToDevice(const unsigned id, // MPI rank
		 const unsigned verbosity,
		 const unsigned geno_count, // # of SNPs
		 const unsigned n, // # of columns in X
		 const unsigned SNPLength, // # length of each SNP
		 double *&d_snptsnp, double *&d_Xtsnp, size_t &d_XtsnpPitch, 
		 double *&d_snpty, char *&d_snpMask, float *&d_f, 
		 double *&d_geno, size_t &d_genoPitch, 
		 double *&d_Xt, size_t &d_XtPitch,
		 double *&d_y,
		 const std::vector<double> &Xty, 
		 const FortranMatrix &XtXi, // G
		 const FortranMatrix &Xt, // only for XtSNP
		 const FortranMatrix &geno,
		 double *&d_nextSNP,
		 const double *y);

void copyUpdateToDevice(unsigned id, unsigned iteration,
			unsigned geno_count, unsigned n,
			char *d_snpMask, 
			int maxFIndex, double *d_Xtsnp, 
			size_t d_XtsnpPitch,
			const FortranMatrix &XtXi,
			const std::vector<double> &Xty);

float getGPUCompTime();

float getGPUMaxTime();

void getMaxFGPU(unsigned id, unsigned iteration, unsigned geno_count, 
	     std::vector<float> &Fval, 
	     unsigned maxFIndex, float *d_f);

// defined in cuda_blas.cu:
void columnDot_gpu(const double *d_mat, unsigned n_rows, uint64_t n_cols, 
		   unsigned columnPitchInWords, // words between tops of columns
		   double *d_result, 
		   unsigned resultStrideInWords // words between result elements
		   );
