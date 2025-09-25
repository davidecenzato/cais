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
Input file smaller than 4.29 GB will take 5n bytes in RAM for the Text and CA plus the overhead
for storing the string boundaries with the sdsl bit vector.
The files larger than 4.29 GB will take 9n bytes plus the overhead for the bitvector.
*/

// algorithms for computing different BWT variants
#include "BWTalgos.h"
// memory counter
#include "external/malloc_count/malloc_count.h"

// function that prints the instructions for using the tool
void print_help(char** argv) {
  std::cout << "Usage: " << argv[ 0 ] << " [options] <input filename>" << std::endl;
  std::cout << "  Options: " << std::endl
        << "\t-b B\tselect one BWT variant to construct (eBWT|dolEBWT|BWT|BBWT) " << std::endl
        << "\t    \teBWT: the original extended BWT of Mantaci et al. " << std::endl
        << "\t    \tdolEBWT: the eBWT where each string is terminated with a dollar " << std::endl
        << "\t    \tBWT: the BWT of a string without the terminating dollar " << std::endl
        << "\t    \tBBWT: the original Bijective BWT of Gil and Scott " << std::endl
        << "\t-o O\tbasename for the output files, def. <output filename>" << std::endl
        << "\t-i I\tselect the input format (FASTA|FASTQ|plain_text) " << std::endl
        << "\t-p  \trun the algorithm variant for datasets of primitive strings, def. False " << std::endl 
        << "\t-u  \tuse an uncompressed bitvector to represent string boundaries, def. False " << std::endl 
        << "\t-c  \twrite the conjugate array, def. False " << std::endl
        << "\t-d  \twrite the document array (only for eBWT and dolEBWT), def. False " << std::endl 
        << "\t-v  \tset verbose mode, def. False " << std::endl;
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
  while ((c = getopt( argc, argv, "db:i:pcuo:vh") ) != -1) {
    switch(c) {
      case 'b':
        sarg.assign( optarg );
        arg.variant.assign( sarg ); break;
        // set the BWT variant to compute
      case 'i':
        sarg.assign( optarg );
        arg.format.assign( sarg ); break;
        // set the input format
      case 'o':
        sarg.assign( optarg );
        arg.outname.assign( sarg ); break;
        // store the output files path
      case 'p':
        arg.no_powers = true; break;
        // algorithm variant for inputs wo. powers
      case 'u':
        arg.sparse = false; break;
        // use sparse bitvector
      case 'c':
        arg.ca = true; break;
        // write the conjugate array
      case 'd':
        arg.doc = true; break;
        // write the document array
      case 'v':
        arg.verbose = true; break;
        // verbose mode
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
  if(arg.variant == ""){
    std::cerr << "Error, select a BWT variant." << std::endl;
    exit(1);
  }
}

int main(int argc, char** argv)
{
  // translate command line arguments
  Args arg;
  parseArgs(argc, argv, arg);
  // select which BWT variant to compute
  if(arg.variant == "eBWT"){
    if(arg.verbose) std::cout << "Computing the eBWT of: " << arg.filename << std::endl;
    if(arg.sparse){ compute_ebwt_sparse(arg, false); }
    else{ compute_ebwt_plain(arg, false); }
  }
  else if(arg.variant == "dolEBWT"){
    if(arg.verbose) std::cout << "Computing the dollar eBWT of: " << arg.filename << std::endl;
    if(arg.sparse){ compute_ebwt_sparse(arg, true); }
    else{ compute_ebwt_plain(arg, true); }
  }
  else if(arg.variant == "BWT"){
    if(arg.verbose) std::cout << "Computing the BWT of: " << arg.filename << std::endl;
    compute_bwt_wo_dol(arg);
  }
  else if(arg.variant == "BBWT"){
    if(arg.verbose) std::cout << "Computing the BBWT of: " << arg.filename << std::endl;
    if(arg.sparse){ compute_bbwt_sparse(arg); }
    else{ compute_bbwt_plain(arg); }
  }
  else{
    std::cout << "Error, select a valid BWT variant... exiting." << std::endl;
    exit(1);
  }

  return 0;
}