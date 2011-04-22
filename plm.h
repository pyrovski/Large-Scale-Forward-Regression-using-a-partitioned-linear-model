#include "type.h"
#include "fortran_matrix.h"
#define iterationLimit 50
unsigned plm_GPU(unsigned, unsigned, unsigned, unsigned, ftype*, ftype*, unsigned, ftype,
	unsigned, ftype*, unsigned*, ftype*);
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
		 ftype *&d_snptsnp, ftype *&d_Xtsnp, unsigned &d_XtsnpPitch, 
		 ftype *&d_snpty, unsigned *&d_snpMask, ftype *&d_f,
		 const vector<double> &SNPtSNP, const FortranMatrix &XtSNP,
		 const vector<double> &SNPty,
		 const vector<double> &Xty, const FortranMatrix &XtXi, 
		 const vector<unsigned> &snpMask);
void copyFromDevice();
void copyUpdateToDevice();
float getGPUCompTime();
float getGPUMaxTime();
void getMaxF();
