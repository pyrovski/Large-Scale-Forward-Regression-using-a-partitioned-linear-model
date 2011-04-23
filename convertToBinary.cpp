#include <iostream>
#include <fstream>

#include <unistd.h>

using namespace std;

const char usage1[] = "usage: ",
	  usage2[] = " <input file> <output file>";

int
main(int argc, char *argv[]){
  if(argc < 3)
    cout << usage1 << argv[0] << usage2;
}
