# cais
cais is a tool to compute four BWT variants using the conjugate array induced sorting algorithm (cais).

# Usage

```
Usage: ./cais [options] <input filename>
  Options: 
  -b B  select one BWT variant to construct (eBWT|dolEBWT|BWT|BBWT) 
        eBWT:    the original extended BWT of Mantaci et al. 
        dolEBWT: the eBWT where each string is terminated with a dollar 
        BWT:     the BWT of a string without the terminating dollar 
        BBWT:    the original Bijective BWT of Gil and Scott 
  -o O  basename for the output files, def. <input filename>
  -i I  select the input format (FASTA|FASTQ|plain_text) 
  -p    run the algorithm variant for datasets of primitive strings only, def. False 
  -u    use an uncompressed bitvector to represent string boundaries, def. False 
  -c    write the conjugate array, def. False 
  -d    write the document array (only for eBWT and dolEBWT), def. False 
  -v    set verbose mode, def. False
```
When computing the eBWT and the dolEBWT, you can choose the format of your input between FASTA and FASTQ.
In addition, with the -u flag, the tool uses an uncompressed bit vector to store the strings' boundaries. Without the -p flag, the eBWT can be computed for any string collection, including power words (strings of the form u^k with k > 1); with the -p flag, the tool runs a simpler algorithm variant that assumes all strings are primitive, and with power words in the input its output is not correct.
When computing the BWT and the BBWT, the input is taken as one single text. The input cannot contain characters with byte value 0 or 1, which are used internally by the algorithm.

### Requirements

The cais tool requires:
* A modern C++14 compiler, such as `g++` version 4.9 or higher.
* A 64-bit macOS or Linux system.

# Example

### Download and Compile

```console
git clone https://github.com/davidecenzato/cais.git
cd cais
git submodule update --init --recursive

mkdir external/sdsl-lite/installed
external/sdsl-lite/install.sh external/sdsl-lite/installed/

make
```

### Run on Example Data

```console
# Construct the eBWT, the GCA and the document array of a toy data set containing power words
./cais -b eBWT -i FASTA -c -d data/powers.fasta
```

### Output files format

The `./cais` executable generates different output files depending on the input flag configuration. You can modify the output file prefix using the `-o` flag. If your input files are larger than 4GB you are required to use the `./cais64` executable.

- **BWT files**:  
  Eight files containing one of the four BWT variants computed by this software, along with the indexes necessary to invert the output transforms. The BWT files are stored in ASCII format (1 byte per character), while the index files are stored in binary format (4 or 8 bytes per string, using `./cais` and `./cais64`, respectively):  
  - `.ebwt`, `.ei` and `.exp`: the extended BWT of Mantaci et al. (using `-b eBWT`). For each input string, `.ei` contains the position in the eBWT of its first character, and `.exp` its exponent (`.exp` is not written with `-p`). To invert the eBWT, follow the LF-mapping from each position in `.ei` until it comes back to that position: this gives the primitive root of the string, which is then repeated `.exp` times.
  - `.dolebwt` and `.di`: the extended BWT of a string collection where all sequences are terminated with a dollar character (using `-b dolEBWT`).  
  - `.bwt` and `.i`: the BWT of a single text without the final terminating dollar (using `-b BWT`). 
  - `.bbwt` and `.bi`: the Bijective BWT of Gil and Scott (using `-b BBWT`).

- **Additional**:  
  Two files storing the (Generalized) Conjugate Array and the Document array in binary format (4 or 8 bytes per entry, using `./cais` and `./cais64`, respectively): 
  - `.gca`: the Generalized Conjugate Array of the input collection (using `-c`). When the `-d` flag is not set, the indexes are given w.r.t. the concatenation of the input strings. Conjugates of power words are ordered by the exponent of their root, then by input order.
  - `.ca`: the Conjugate Array of one input text (using `-c`).
  - `.da`: the Document Array for the GCA (using `-d` and `-c` together).

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

If you notice any bugs, please feel free to report them by opening a Git issue or by contacting us at davide[x]cenzato[y]unive[x]it, where x = '.' and y = '@'.