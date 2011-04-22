#include "type.h"
int plm_GPU(unsigned, unsigned, unsigned, unsigned, ftype*, ftype*, unsigned, ftype,
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

int copyToDevice();
int copyFromDevice();
int copyUpdateToDevice();
float getGPUCompTime();
float getGPUMaxTime();
int getMaxF();
