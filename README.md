# cais
cais is a tool that computes four BWT variants using the conjugate array induced sorting algorithm (cais).

# Usage

```
Usage: ./cais [options] <input filename>
  Options:
        -e      construct the extended BWT (eBWT), def. True
        -d      construct the dollar eBWT (dolEBWT), def. False
        -b      construct the BWT of the text without dollar (BWT), def. False 
        -t      construct the bijective BWT (BBWT), def. False
        -f      take in input a fasta file (only for eBWT and dolEBWT), def. True
        -q      take in input a fastq file (only for eBWT and dolEBWT), def. False
        -s      use sparse bitvector (less space with long strings), def. False
        -c      write the conjugate array, def. False
        -a      write the document array (only for eBWT and dolEBWT), def. False
        -v      set verbose mode, def. False
        -o O    basename for the output files, def. <input filename>
```
When computing the eBWT and eBWT you can choose the format of your input between fasta and fastq format.
When computing the BWT and the BBWT the input is taken as a single text.
With the -s flag, the tool uses a sparse bit vector to store the strings' boundaries, with the -a flag the tool
writes the generalized conjugate array in two files with extensions '.da' and '.gca' containing the string identifier and
the offset of each conjugate.


### Requirements

The cais tool requires:
* A modern C++11 compiler such as `g++` version 4.9 or higher.

# Example

### Download and Compile

```console
git clone https://github.com/davidecenzato/cais.git
cd cais
git submodule update --init --recursive

cd external/sdsl-lite
mkdir installed
./install.sh installed/

cd ../../
make
```

### Run on Example Data

```console
// Construct the eBWT and the gCA on a toy data set
./cais test.fasta -e -c -f 
```

### Output files format

The `./cais` executable generates different output files depending on the input parameter configuration. You can modify the output file prefix using the `-o` flag. If your input files are larger than 4GB you can use the `./cais64` executable.

- **BWT variant files**:  
  Eight files containing one of the four BWT variants computed by this software, along with the indexes necessary to invert the computed transforms. The BWT files are stored in ASCII format (1 byte per character), while the index files are stored in binary format (4 and 8 bytes per string, using `./cais` and `./cais64`, respectively):  
  - `.ebwt` and `.ei`: the extended BWT of Mantaci et al. (using `-e` flag). The `.ei` file can be used to invert the eBWT iff. the input strings are primitive.   
  - `.dolebwt` and `.di`: the extended BWT of a collection of strings where all sequences are terminated with a dollar character (using `-d` flag).  
  - `.bwt` and `.i`: the BWT of a single text without the final dollar (using `-b` flag). 
  - `.bbwt` and `.bi`: the Bijective BWT of Gill and Scott (using `-t` flag).  

- **Additional**:  
  Two files storing the (Generalized) Conjugate Array and the Document array (4 and 8 bytes per character, using `./cais` and `./cais64`, respectively): 
  - `.gca`: the Generalized Conjugate Array of the input collection (using `-c` togheter with `-e`, `-d` and `-t` flags). When the `-a` flag is not used, the indexes are given w.r.t. the concatenation of the input strings.
  - `.ca`: the Conjugate Array of the input text (using `-c` togheter with `-b` flag).
  - `.da`: the Document Array for the GCA (using `-c` and `-a` togheter with `-e` and `-d` flags).

# External resources

* [malloc_count](https://github.com/bingmann/malloc_count)
* [sdsl-lite](https://github.com/simongog/sdsl-lite)

# Citation 

[1] Christina Boucher, Davide Cenzato, Zsuzsanna Lipták, Massimiliano Rossi, Marinella Sciortino: *Computing the Original eBWT Faster, Simpler, and with Less Memory*. 28th International Symposium on String Processing and Information Retrieval (SPIRE 2021). Lecture Notes in Computer Science vol. 12944: 129-142 (2021) ([go to the paper](https://link.springer.com/chapter/10.1007/978-3-030-86692-1_11)).

[2] Christina Boucher, Davide Cenzato, Zsuzsanna Lipták, Massimiliano Rossi, Marinella Sciortino:
*Computing the original eBWT faster, simpler, and with less memory*. CoRR abs/2106.11191 (2021). ([go to the paper](https://arxiv.org/abs/2106.11191))

If you use this tool in an academic setting, please cite this work as follows:

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

### Contacts:

If you notice any bugs, please feel free to report them by opening a Git issue or by contacting us at davidecenzato Unive email.