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
