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
using namespace std;

const char usage1[] = "usage: ",
	  usage2[] = " <input file> <output file>";

int
main(int argc, char *argv[]){
  if(argc < 3)
    cout << usage1 << argv[0] << usage2;


  FILE *outfile = fopen(argv[2], "w");
  if(!outfile){
    cerr << "failed to open output file" << endl;
    exit(1);
  }
  long flags = fcntl(fileno(outfile), F_GETFL);
  if(fcntl(fileno(outfile), F_SETFL, flags | O_NONBLOCK)){
    perror("fcntl");
    exit(1);
  }
  ifstream instream(argv[1]);
  const unsigned arrayLength = 1024*1024;
  double val[arrayLength];
  uint64_t count = 0;
  while(!instream.eof() && !instream.fail()){
    unsigned i;
    for(i = 0; i < arrayLength; i++){
      //int status = fscanf(infile, "%le", &val[i]);
      //if(status != 1)
      //break;
      instream >> val[i];
      count++;
      if(instream.eof() || instream.fail())
	break;
    }
    
    fwrite(val, sizeof(double), i, outfile);
  }
  if(instream.fail() && !instream.eof()){
    cerr << "file error on input file" << endl;
    cerr << "read " << count << " entries" << endl;
  }
  fclose(outfile);
  instream.close();
  return 0;
  
}
