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
using namespace std;

const char usage1[] = "usage: ",
	  usage2[] = "-r <input rows> -c <input cols> <input file> -o <output file>";

void usage(char ** argv){
  cout << usage1 << argv[0] << usage2;
}

const int blockRows = 128;
const int blockCols = 128;

int
main(int argc, char *argv[]){
  if(argc < 5){
    usage(argv);
    exit(1);
  }
  
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
      usage(argv);
      exit(1);
    }
  }

  if(!rows || !cols || !inputFilename){
    usage(argv);
    exit(1);
  }

  

  FILE *inFile = 0;

  if(inPlace){
    inFile = fopen(inputFilename, "w+");
  }
  
  assert(inFile);
  int ifd = fileno(inFile);
  assert(ifd > 0);
  
  struct stat fileStat;
  status = fstat(ifd, &fileStat);
  assert(!status);
  
  double *inData = 0,
    *outData = 0;
  
  if(inPlace){
    inData = (double*)
      mmap(0, fileStat.st_size, PROT_READ | PROT_WRITE, MAP_SHARED,
	   ifd, 0);
    assert(inData);
    outData = inData;
  }

  for(int i = 0; i < rows; i += blockRows){
    int kIter = i + blockRows <= rows ? blockRows : rows - i;
    for(int j = 0; j < cols; j += blockCols){
      int lIter = j + blockCols <= cols ? blockCols : cols - j;
      for(int k = 0; k < kIter; k++)
	for(int l = 0; l < lIter; l++)
	  outData[j * cols + i] = inData[i * cols + j];
    }
  }

  if(inPlace)
    status = munmap(inData, fileStat.st_size);
  
  return 0;
  
}
