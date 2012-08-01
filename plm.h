/*


Copyright (c) 2011, The Arizona Board of Regents on behalf of 
The University of Arizona

All rights reserved.

Developed by Peter Bailey, Tapasya Patki, and Greg Striemer with
support from the iPlant Collaborative as a collaboration between
participants at BIO5 at The University of Arizona (the primary hosting
institution), Cold Spring Harbor Laboratory, The University of Texas
at Austin, and individual contributors. Find out more at
http://www.iplantcollaborative.org/.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are
met:

 * Redistributions of source code must retain the above copyright 
   notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright 
   notice, this list of conditions and the following disclaimer in the 
   documentation and/or other materials provided with the distribution.
 * Neither the name of the iPlant Collaborative, BIO5, The University 
   of Arizona, Cold Spring Harbor Laboratory, The University of Texas at 
   Austin, nor the names of other contributors may be used to endorse or 
   promote products derived from this software without specific prior 
   written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


*/
#include <vector>

#include "type.h"
#include "fortran_matrix.h"

/*! @todo make hardIterationLimit a hard limit, and iterationLimit a soft limit.
  hardIterationLimit is constrained by constant memory size, which varies
  according to device family, so we can define hardIterationLimit as 
  a macro function of 
 */
//#define fixedPlusIteration_limit 89

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
