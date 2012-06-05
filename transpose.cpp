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
#include <iostream>
#include <fstream>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>
#include <algorithm>
using namespace std;

const char usage1[] = "usage: ",
	  usage2[] = "-r <input rows> -c <input cols> <input file> -o <output file> [-i (for square matrices only!)]\nAssuming input is row-major\n";

void usage(char ** argv){
  cerr << usage1 << argv[0] << " " << usage2;
}

#ifndef blockRows
#define blockRows 8192
#endif

#ifndef blockCols
#define blockCols blockRows
#endif

#define min(a,b) ((a) < (b) ? (a) : (b))

int
main(int argc, char *argv[]){
 
  int status;
  int rows = 0, cols = 0;
  char *outputFilename = 0, *inputFilename = 0;
  bool inPlace = false;

  while((status = getopt(argc, argv, "r:c:o:i")) != -1){
    switch(status){
    case 'r':
      rows = strtoul(optarg, 0, 0);
      break;
    case 'c':
      cols = strtoul(optarg, 0, 0);
      break;
    case 'o':
      outputFilename = optarg;
      break;
    case 'i':
      inPlace = true;
      break;
    default:
      cerr << "unexpected parameter: " << status << endl;
      usage(argv);
      exit(1);
    }
  }

  
  if(optind >= argc){
    cerr << "must provide input file" << endl;
    usage(argv);
    exit(1);
  }

  if(inPlace)
    cols = rows;
  
  inputFilename = argv[optind];

  if(!rows || !cols){
    cerr << "must provide rows, columns" << endl;
    usage(argv);
    exit(1);
  }

  FILE *inFile = 0, *outFile = 0;
  int ifd = -1, ofd = -1;
  struct stat fileStat;
  double *inData = 0,
    *outData = 0;

  if(inPlace)
    inFile = fopen(inputFilename, "r+");
  else
    inFile = fopen(inputFilename, "r");

  assert(inFile);
  ifd = fileno(inFile);
  assert(ifd > 0);
  
  status = fstat(ifd, &fileStat);
  assert(!status);
    
#ifdef _DEBUG
  cerr << "size: " << fileStat.st_size << endl;
#endif

  if(inPlace){
    inData = (double*)
      mmap(0, fileStat.st_size, PROT_READ | PROT_WRITE, 
	   MAP_SHARED | MAP_NORESERVE,
	   ifd, 0);
  } else {
    inData = (double*)
      mmap(0, fileStat.st_size, PROT_READ, 
	   MAP_SHARED | MAP_NORESERVE,
	   ifd, 0);
  }
    
  if(inData == (void*)-1){
    //! @todo if mmap fails, do it out of core
    perror("mmap input failed");
    exit(2);
  }

  if(inPlace){
    outData = inData;
  } else {
    outFile = fopen(outputFilename, "w+");
    ofd = fileno(outFile);
    assert(ofd > 0);

    status = ftruncate(ofd, fileStat.st_size);
    assert(!status);
    
    outData = (double*)
      mmap(0, fileStat.st_size, PROT_WRITE, 
	   MAP_SHARED | MAP_NORESERVE | MAP_POPULATE,
	   ofd, 0);
  
    if(outData == (void*)-1){
      //! @todo if mmap fails, do it out of core
      perror("mmap output failed");
      exit(2);
    }
  } 
  
  if(inPlace){
    double *swapData1 = (double*)malloc(blockRows * blockCols * sizeof(double));
    double *swapData2 = (double*)malloc(blockRows * blockCols * sizeof(double));
    
    for(int i = 0; i < rows; i += blockRows){
      int kIter = i + blockRows <= rows ? blockRows : rows - i;
      int jLim = min(i, cols);
      for(int j = 0; j < jLim; j += blockCols){
	int lIter = j + blockCols <= jLim ? blockCols : jLim - j;
	for(int k = 0; k < kIter; k++){
	  for(int l = 0; l < lIter; l++){
	    swapData1[l * blockCols + k] = 
	      inData[(i + k) * cols + (j + l)];
	    // where will this block end up?
	    swapData2[k * blockCols + l] = 
	      inData[(j + l) * rows + (i + k)];
	  }
	}
	/* Copy transposed data 1 to output.
	   j + l = rows of output
	   i + k = cols of output
	*/
	for(int l = 0; l < lIter; l++)
	  memcpy(outData + (j + l) * rows + i, 
		 swapData1 + l * blockCols, 
		 kIter * sizeof(double));

	/* Copy transposed data 2 to output.
	   j + l = rows of output
	   i + k = cols of output
	*/
	for(int k = 0; k < kIter; k++)
	  memcpy(outData + (i + k) * cols + j, 
		 swapData2 + k * blockCols, 
		 lIter * sizeof(double));
      }
    }
    // handle blocks on the diagonal
    int iLim = min(rows, cols);
    for(int i = 0; i < iLim; i += blockRows){
      int kIter = i + blockRows <= rows ? blockRows : rows - i;
      int lIter = i + blockCols <= cols ? blockCols : cols - i;
      for(int k = 0; k < kIter; k++)
	for(int l = k + 1; l < lIter; l++)
	  swap(inData[(i + k) * cols + (i + l)], 
	       outData[(i + l) * rows + i + k]);
    }
    free(swapData1);
    free(swapData2);
  } else {
    for(int i = 0; i < rows; i += blockRows){
      int kIter = i + blockRows <= rows ? blockRows : rows - i;
      for(int j = 0; j < cols; j += blockCols){
	int lIter = j + blockCols <= cols ? blockCols : cols - j;
	for(int k = 0; k < kIter; k++)
	  for(int l = 0; l < lIter; l++)
	    outData[(j + l) * rows + i + k] = 
	      inData[(i + k) * cols + (j + l)];
      }
    }
  }

  status = munmap(inData, fileStat.st_size);
  if(!inPlace){
    status = munmap(outData, fileStat.st_size);
    fclose(outFile);
  }

  fclose(inFile);
  return 0;
  
}
