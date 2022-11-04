# cais
cais is a tool that computes several BWT variants using the conjugate array induced sorting algorithm (cais).

# Usage

```
Usage: ./cais <input filename> [options]
  Options:
        -e      construct the extended BWT (eBWT), def. True
        -d      construct the dollar eBWT (dolEBWT), def. False
        -b      construct the BWT without dollar (BWT), def. False
        -t      construct the bijective BWT (BBWT), def. False
        -f      take in input a fasta file (only for eBWT and dolEBWT), def. True
        -q      take in input a fastq file (only for eBWT and dolEBWT), def. False
        -s      write the conjugate array, def. False
        -v      set verbose mode, def. False
        -o O    basename for the output files, def. <input filename>

```
When computing the eBWT and eBWT you can choose the format of your input between fasta and fastq format.
When computing the BWT and the BBWT the input is taken as a single text.

### Requirements

The cais tool requires:
* A modern C++11 compiler such as `g++` version 4.9 or higher.
* The sdsl-lite library installed.

# Example

### Download and Compile

```console
https://github.com/davidecenzato/cais.git
cd cais
git submodule update --init --recursive
make
```

### Run on Example Data

```console
// Construct the eBWT and the gCA on a toy data set
./cais test.fasta -e -s -f 
```
# External resources

* [malloc_count](https://github.com/bingmann/malloc_count)
* [sdsl-lite](https://github.com/simongog/sdsl-lite)

# Citation 

If you use this tool in an academic setting, please cite the following paper, in which the algorithm was developed for computing of the eBWT and the BWT:

[1] Christina Boucher, Davide Cenzato, Zsuzsanna Lipták, Massimiliano Rossi, Marinella Sciortino: *Computing the Original eBWT Faster, Simpler, and with Less Memory*. 28th International Symposium on String Processing and Information Retrieval (SPIRE 2021). Lecture Notes in Computer Science vol. 12944: 129-142 (2021).

### cais
    @inproceedings{BoucherCL0S21a,
      author    = {Christina Boucher and
                   Davide Cenzato and
                   {\relax Zs}uzsanna Lipt{\'{a}}k and
                   Massimiliano Rossi and
                   Marinella Sciortino},
      title     = {Computing the Original e{BWT} Faster, Simpler, and with Less Memory},
      booktitle = {Proceedings of 28th International Symposium in String Processing and Information Retrieval {SPIRE} 2021},
      series    = {Lecture Notes in Computer Science},
      volume    = {12944},
      pages     = {129--142},
      year      = {2021}
    }

# Authors

### Theoretical results:

* Christina Boucher
* Davide Cenzato
* Zsuzsanna Lipták
* Massimiliano Rossi
* Marinella Sciortino

### Implementation:

* [Davide Cenzato](https://github.com/davidecenzato) 