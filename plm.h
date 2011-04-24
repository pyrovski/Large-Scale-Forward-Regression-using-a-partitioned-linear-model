#include "type.h"
#include "fortran_matrix.h"
#define iterationLimit 50
unsigned plm_GPU(unsigned geno_count, unsigned blockSize, unsigned sharedSize, 
		 unsigned m, ftype* d_snptsnp, ftype* d_Xtsnp, 
		 unsigned d_XtsnpPitch, ftype ErrorSS, unsigned V2, 
		 ftype* d_snpty, unsigned* d_snpMask, ftype* d_f) throw(int);
/*
geno_count, n, n * sizeof(ftype), 
	m ,        
	d_snptsnp, 
	d_Xtsnp, 
	d_XtsnpPitch, 
	glm_data.ErrorSS, glm_data.V2, 
	d_snpty, 
	d_snpMask,
	d_f
*/

void copyToDevice(const unsigned geno_count, const unsigned n, 
		 ftype *&d_snptsnp, ftype *&d_Xtsnp, size_t &d_XtsnpPitch, 
		 ftype *&d_snpty, unsigned *&d_snpMask, ftype *&d_f,
		  const std::vector<double> &SNPtSNP, const FortranMatrix &XtSNP,
		 const std::vector<double> &SNPty,
		 const std::vector<double> &Xty, const FortranMatrix &XtXi, 
		 const std::vector<unsigned> &snpMask);
void copyUpdateToDevice(unsigned geno_count, unsigned n,
		       unsigned *d_snpMask, 
		       unsigned maxFIndex, ftype *d_Xtsnp, 
		       size_t d_XtsnpPitch,
			const std::vector<unsigned> &snpMask,
		       FortranMatrix &XtSNP, const FortranMatrix &XtXi,
		       const std::vector<double> &Xty);
float getGPUCompTime();
float getGPUMaxTime();
void getMaxF(unsigned id, unsigned iteration, unsigned geno_count, 
	     std::vector<double> &Fval, 
	     unsigned maxFIndex, ftype *d_f);
