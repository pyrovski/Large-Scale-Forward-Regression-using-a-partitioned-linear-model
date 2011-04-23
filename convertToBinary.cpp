#include <iostream>
#include <fstream>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <stdlib.h>

using namespace std;

const char usage1[] = "usage: ",
	  usage2[] = " <input file> <output file>";

int
main(int argc, char *argv[]){
  if(argc < 3)
    cout << usage1 << argv[0] << usage2;

  FILE *infile = fopen(argv[1], "r");
  FILE *outfile = fopen(argv[2], "w");
  if(!infile || !outfile){
    cerr << "failed to open input or output file" << endl;
    exit(1);
  }
  long flags = fcntl(fileno(outfile), F_GETFL);
  if(fcntl(fileno(outfile), F_SETFL, flags | O_NONBLOCK)){
    perror("fcntl");
    exit(1);
  }
  while(!feof(infile) && !ferror(infile)){
    double val;
    int status = fscanf(infile, "%le", &val);
    if(status == 1)
      fwrite(&val, sizeof(double), 1, outfile);
  }
  fclose(outfile);
  fclose(infile);
  return 0;
  
}
