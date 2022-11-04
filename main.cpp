#include <string>
#include <iostream>
#include <getopt.h>

/*
Description:
Tool for computing different BWT variants: extended BWT (eBWT), dollar eBWT (dolEBWT),
the BWT without dollar (BWT) and the Bijective BWT (BBWT) using the Conjugate array 
induced sorting (cais) algorithm.
It can also give in output the generalized conjugate array of a string collection
or the conjugate array of a text
   "Computing the Original eBWT Faster, Simpler, and with Less Memory"
   by Christina Boucher, Davide Cenzato, Zsuzsanna Lipták, Massimiliano Rossi and Marinella Sciortino
   Proc. of SPIRE '21, doi = {10.1007/978-3-030-86692-1\_11}
   
The input file cannot contain the characters < 2 which are used internally by the algorithm.
Input file smaller than 4.29 GB will take 5n + n/8 bytes in RAM, while input files larger than
4.29 GB will take 9n + n/8 bytes.
*/

// algorithms for computing different BWT variants
#include "BWTalgos.h"
// memory counter
#include "external/malloc_count/malloc_count.h"

// function that prints the instructions for using the tool
void print_help(char** argv) {
  std::cout << "Usage: " << argv[ 0 ] << " <input filename> [options]" << std::endl;
  std::cout << "  Options: " << std::endl
        << "\t-e \tconstruct the extended BWT (eBWT), def. True " << std::endl
        << "\t-d \tconstruct the dollar eBWT (dolEBWT), def. False " << std::endl
        << "\t-b \tconstruct the BWT without dollar (BWT), def. False " << std::endl
        << "\t-t \tconstruct the bijective BWT (BBWT), def. False " << std::endl
        << "\t-f \ttake in input a fasta file (only for eBWT and dolEBWT), def. True " << std::endl
        << "\t-q \ttake in input a fastq file (only for eBWT and dolEBWT), def. False " << std::endl 
        //<< "\t-n \ttake in input a new line sep file, def. False " << std::endl to be implemented 
        << "\t-s \twrite the conjugate array, def. False " << std::endl
        << "\t-v \tser verbose mode, def. False " << std::endl
        << "\t-o O\tbasename for the output files, def. <input filename>" << std::endl;

  exit(-1);
}
// function for parsing the input arguments
void parseArgs( int argc, char** argv, Args& arg ) {
  int c;
  extern int optind;

  puts("==== Command line:");
  for(int i=0;i<argc;i++)
    printf(" %s",argv[i]);
  puts("");
 
  std::string sarg;
  while ((c = getopt( argc, argv, "edbtfqso:vh") ) != -1) {
    switch(c) {
      case 'e':
        arg.variant = 0; break;
        // compute the eBWT
      case 'd':
        arg.variant = 1; break;
        // compute the dolEBWT
      case 'b':
        arg.variant = 2; break;
        // compute the BWT
      case 't':
        arg.variant = 3; break;
        // compute the BBWT
      case 'f':
        arg.format = 1; break;
        // take in input a fasta file
      case 'q':
        arg.format = 2; break;
        // take in input a fastq file
      case 'v':
        arg.verbose = true; break;
        // verbose mode
      case 's':
        arg.ca = true; break;
        // write the suffix array
      case 'o':
        sarg.assign( optarg );
        arg.outname.assign( sarg ); break;
        // store the output files path
      case 'h':
        print_help(argv); exit(-1);
        // fall through
      default:
        std::cout << "Unknown option. Use -h for help." << std::endl;
        exit(-1);
    }
  }
  // the only input parameter is the file name
  if (argc == optind+1) {
    arg.filename.assign( argv[optind] );
  }
  else {
    std::cout << "Invalid number of arguments" << std::endl;
    print_help(argv);
  }
  // set output files basename
  if(arg.outname == "") arg.outname = arg.filename;
  // check BWT variant
  if(arg.variant == -1){
    std::cerr << "Error, select a BWT variant." << std::endl;
    exit(1);
  }
  // check input type
  if(arg.variant > 1){
    // force input format
    if(arg.format > 0) arg.format = 0;
  }
  else{
    //force input format
    if(arg.format == 0) arg.format = 1;
  }
}

int main(int argc, char** argv)
{
  // translate command line arguments
  Args arg;
  parseArgs(argc, argv, arg);
  // select which BWT variant to compute
  switch(arg.variant) {
    case 0:
      if(arg.verbose) std::cout << "Computing the eBWT of: " << arg.filename << std::endl;
      compute_ebwt(arg, false);
      break;
    case 1:
      if(arg.verbose) std::cout << "Computing the dollar eBWT of: " << arg.filename << std::endl;
      compute_ebwt(arg, true);
      break;
    case 2:
      if(arg.verbose) std::cout << "Computing the BWT of: " << arg.filename << std::endl;
      compute_bwt_wo_dol(arg);
      break;
    case 3:
      if(arg.verbose) std::cout << "Computing the BBWT of: " << arg.filename << std::endl;
      compute_bbwt(arg);
      break;
    default:
      std::cout << "Error, select a valid BWT variant... exiting." << std::endl;
      exit(1);
  }

  return 0;
}