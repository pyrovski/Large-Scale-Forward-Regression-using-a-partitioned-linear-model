// System header files
#include <sys/time.h>
#include <string>
#include <fstream>
#include <vector>

// Local project includes
#include "fortran_matrix.h"
#include "glm.h"
#include "GetPot"

using namespace std;

int main()
{
  // Timing variables
  timeval tstart, tstop;

  // Create input file object.  Put the path to your data files here!
  GetPot input_file("reference_glm.in");
  
  // The path on my system to the location of the data files.  Don't forget the trailing
  // slash here, as this will be prepended to the filename below
  string path = "./"; // default path is the current directory.
  path = input_file("path", path.c_str());
  
  // File containing the "population structure".  It is a 4892-by-26 matrix
  string fixed_filename = "fixed.effects.nam.sorted.filtered.dat";

  // File containing the genotypes.  It is a 4892-by-79 matrix.
  string geno_filename = "imputed.marker.chr10.sorted.filtered.dat";

  // File containing the phenotypes.  It is a 4892-by-1 matrix.  The file is
  // arranged in a single column.
  string y_filename = "residuals.chr10.sorted.dat";

  // In Matlab, these sizes are inferred from the data.  In C++, we hard-code them
  // to make reading the data simpler...
  unsigned pop_ind = 4892, fixed_count = 26; // rows, columns of the fixed array
  unsigned geno_ind = 4892, geno_count = 79; // rows, columns of the geno array
  unsigned y_ind = 4892, y_count = 1;        // rows, columns of the y vector

  // Matrix objects for storing the input data
  FortranMatrix fixed(pop_ind, fixed_count);
  FortranMatrix geno(geno_ind, geno_count);
  vector<double> y(y_ind);

  // Begin timing the file IO for all 3 files
  gettimeofday(&tstart, NULL);

  
  // Read the "fixed" array from file
  {
    // Open "fixed" file for reading
    string filename = path + "/" + fixed_filename;
    ifstream fixed_file(filename.c_str());

    if (fixed_file)
      {
	// Loop over all rows and columns, set entries in the fixed matrix
	// double val=99;
	for (unsigned i=0; i<fixed.get_n_rows(); ++i)
	  for (unsigned j=0; j<fixed.get_n_cols(); ++j)
	    {
	      fixed_file >> fixed(i,j);
	      //cout << "val=" << val << endl;
	      //fixed(i,j) = val;
	    }
      }
    else
      {
	cout << "Failed to open file: " << fixed_file << "!!" << endl;
	return 1;
      }
  }

  // Read the geno array from file

  {
    // Open "fixed" file for reading
    ifstream geno_file((path + "/" + geno_filename).c_str());

    if (geno_file)
      {
	// Loop over all rows and columns, set entries in the matrix
	for (unsigned i=0; i<geno.get_n_rows(); ++i)
	  for (unsigned j=0; j<geno.get_n_cols(); ++j)
	    {
	      geno_file >> geno(i,j);
	    }
      }
    else
      {
	cout << "Failed to open file!!" << endl;
	return 1;
      }
  }


  // Read the y-array from file.  Currently stored as a vector since that
  // is how it is passed to the glm function, but could be changed to a
  // FortranMatrix with one column...

  {
    // Open "fixed" file for reading
    ifstream y_file((path + "/" + y_filename).c_str());

    if (y_file)
      {
	// Loop over all rows and columns, set entries in the matrix
	for (unsigned i=0; i<y.size(); ++i)
	  {
	    y_file >> y[i];
	  }
      }
    else
      {
	cout << "Failed to open file!!" << endl;
	return 1;
      }
  }

  gettimeofday(&tstop, NULL);

  {
    // Compute time taken for IO
    const double io_elapsed_time = (static_cast<double>(tstop.tv_sec  - tstart.tv_sec) +
				    static_cast<double>(tstop.tv_usec - tstart.tv_usec)*1.e-6);

    cout << "Time required for I/O: " << io_elapsed_time << " s." << endl;
  }
  
  
  // Debugging
#ifdef _DEBUG
  cout << "fixed effects:" << endl;
  for (unsigned i=0; i<fixed.get_n_rows(); ++i)
    {
      for (unsigned j=0; j<fixed.get_n_cols(); ++j)
	cout << fixed(i,j) << ' ';
      cout << '\n';
    }

  cout << "geno:" << endl;
  for (unsigned i=0; i<geno.get_n_rows(); ++i)
    {
      for (unsigned j=0; j<geno.get_n_cols(); ++j)
	cout << geno(i,j) << ' ';
      cout << '\n';
    }

  cout << "y:" << endl;
  
  for (unsigned i=0; i<y.size(); ++i)
    {
      cout << y[i] << endl;
    }
#endif
  // Version A, Kt a general matrix.
  // Create the Kt matrix.  It has 1 row and fixed_count+2 columns.
  // The entries of Kt are all zero except for the last entry, which is 1.
 FortranMatrix Kt(1,fixed_count+2);
 Kt(0, fixed_count+1) = 1.; // Set last entry = 1

  // Version B, Kt assumed a vector.
//  vector<double> Kt(fixed_count+2);
//  Kt.back() = 1.; // Set last entry = 1


  
  // Debugging:
  // Kt.print("Kt"); 
  
  // An array to hold the results of the GLM calculations
  vector<double> Pval(geno_count);

  // Initialize the X-matrix.  The first column is all ones, the next
  // fixed_count columns are equatl to the fixed matrix, and the last
  // column (which changes) is the i'th column of the geno array.
  FortranMatrix X(geno_ind/*4892*/, fixed_count+2/*28*/);

  // Fill first column of X with 1's
  for (unsigned i=0; i<X.get_n_rows(); ++i)
    X(i,0) = 1.;

  // Fill next fixed_count columns with the fixed array contents
  for (unsigned j=0; j<fixed_count; ++j)
    for (unsigned i=0; i<X.get_n_rows(); ++i)
      X(i,j+1) = fixed(i,j);

  // To hold return values of the GLM call.  Apparently we only use
  // "p" currently.
  GLMData glm_data;

  // Begin timing the computations
  gettimeofday(&tstart, NULL);
  
  // For each column of the geno array, set up the "X" matrix,
  // call the GLM routine, and store the computed p value.  Note
  // we are assuming y is a vector for now (because GLM currently expects
  // a vector for its second argument) but this could be generalized later.
  for (unsigned i=0; i<geno_count; ++i)
    {
      // Fill the last column of X with i'th column of the geno array
      for (unsigned row=0; row<X.get_n_rows(); ++row)
	X(row, fixed_count+1) = geno(row,i);
      
      // Debugging: print X
      // X.print("X");

      // Call the glm function.  Note that X is currently overwritten by this function,
      // and therefore would need to be re-formed completely at each iteration...
      glm(X, y, Kt, glm_data);

      // Debugging: print the computed p value
      // cout << "glm_data.p=" << glm_data.p << endl;

      // Store the computed value in an array
      Pval[i] = glm_data.p;
    }

  // Finish timing the computations
  gettimeofday(&tstop, NULL);

  {
    // Compute time taken for IO
    const double computation_elapsed_time = (static_cast<double>(tstop.tv_sec  - tstart.tv_sec) +
					     static_cast<double>(tstop.tv_usec - tstart.tv_usec)*1.e-6);

    cout << "Time required for computations: "
	      << computation_elapsed_time
	      << " s."
	      << endl;
  }
  
  // Print out the Pval array
cout << "Pval=" << endl;
for (unsigned i=0; i<Pval.size(); ++i)
 cout << Pval[i] << endl;
    
  return 0;
}


// On my Macbook (seconds) timings (algorithm assumes a Kt-matrix). 
// The algorithm which assumes a Kt-vector is only marginally faster, if at all...

// Avg. Time for I/O 
// (1.2631+1.24829+1.29038+1.26401+1.24963)/5 = 1.263082

// Avg. Computation time
// (0.279372+0.282701+0.271749+0.275394+0.274847)/5 = .2768126

// Percentage of total time used for computations:
// .2768126/(1.263082+.2768126) = .17976074
