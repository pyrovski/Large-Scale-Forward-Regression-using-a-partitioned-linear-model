#include "fortran_matrix.h"

int pinv(FortranMatrix &M, std::vector<double> &rhs, FortranMatrix &Mi, 
	 std::vector<double> &beta, double &tol, int id);
