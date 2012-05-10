#include <iostream>
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

int main(int argc, char **argv){
  if(argc < 2){
    cerr << "must provide input file" << endl;
    exit(1);
  }
  
  FILE * file = fopen(argv[1], "r");

  assert(file);

  uint64_t elem = 0, t_elem = 0;

  const int bufLength = 1024*1024;
  double buf[bufLength];

  while(!feof(file) && !ferror(file)){
    t_elem = fread(buf, sizeof(double), bufLength, file);
    for(int i = 0; i < t_elem; i++)
      if(buf[i] < 0.0){
	cerr << "negative element at index " 
	     << elem + i << ": " << buf[i] << endl;
	exit(1);
      }
    elem += t_elem;
  }
  cout << "read " << elem << " elements (" 
       << elem * sizeof(double) << " bytes)" << endl;
}
