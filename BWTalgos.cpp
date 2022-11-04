#include "BWTalgos.h"

// function for constructing and storing the eBWT and the generalized
// conjugate array of a string collection
void compute_ebwt(Args arg, bool concat){ 
    // input text initialization
    std::vector<unsigned char> Text;
    uint_s n=0, ns=0;
    std::vector<uint_s> onset = {0};
    // initialize gap encoded bit-vector 
    sdsl::sd_vector<> b;
    sdsl::sd_vector<>::rank_1_type r_s;
    sdsl::sd_vector<>::select_1_type s_s;
    // initalize output files
    std::string ebwtfile, ifile, cafile;
    FILE *ebwt, *I, *gCa;
    // initialize output basename
    std::string basename = arg.outname;
    // initialize boolean values
    int format = arg.format;
    bool printca = arg.ca;
    // read input
    if(format==1){
        if(arg.verbose) std::cout << "load fasta file" << std::endl;
        // load a fasta file
        load_fasta(arg.filename.c_str(),Text,onset,n,ns,concat);
    }else{
        if(arg.verbose) std::cout << "load fastq file" << std::endl;
        // load a fastq file
        load_fastq(arg.filename.c_str(),Text,onset,n,ns,concat);
    }
    // construct gap-encoded bit vector marking the starting position
    // of each string in Text
    sdsl::sd_vector_builder builder(n+1,onset.size());
    for(auto idx: onset){ builder.set(idx); }
    b = sdsl::sd_vector<>(builder);
    // clear onset vector
    onset.clear();
    // initialize the circular Suffix array
    std::vector<uint_s> CA(n,0);
    // initialize the names of the output files
    // compute circular suffix array
    auto start = std::chrono::steady_clock::now();
    cais(&Text[0], &CA[0], n, alph_size, b);
    auto end = std::chrono::steady_clock::now();
    if(arg.verbose){
        std::cout << "gCA construction: Elapsed time in seconds: "
        << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
        << " s" << std::endl;
    }
    start = std::chrono::steady_clock::now();
    // initialize rank and select data structures
    r_s = sdsl::sd_vector<>::rank_1_type(&b);
    s_s = sdsl::sd_vector<>::select_1_type(&b);
    // initialize the names of the output files
    if(concat){
        // select the dollarebwt extensions
        ebwtfile = basename + std::string(".dolebwt");
        ifile = basename + std::string(".di");        
    }else{
        // select the ebwt extensions
        ebwtfile = basename + std::string(".ebwt");
        ifile = basename + std::string(".i");
    }
    // open output files
    if((ebwt = fopen(ebwtfile.c_str(), "w")) == nullptr){
        std::cerr << "open() file " + ebwtfile + " failed!" << std::endl;
        exit(-1);
    }
    if((I = fopen(ifile.c_str(), "w")) == nullptr){
        std::cerr << "open() file " + ifile + " failed!" << std::endl;
        exit(-1);
    }
    // write the ebwt and the i vector to file
    for(size_t i=0; i<CA.size(); ++i){
        char c;
        // if cSA[i] is not the first position of a string
        if(b[CA[i]]!=1){ c = Text[CA[i]-1];  }
        // write the position in the i vector
        else{ 
            c = Text[s_s(r_s(CA[i]+1)+1)-1];
            if(fwrite(&i,sizeof(uint32_t),1,I)!=1){
                std::cerr << "I file write error... exiting!" << std::endl;
                exit(-1);
            }
          }
        // write a new character in the ebwt
        putc(c, ebwt);
    }
    // close output files
    fclose(ebwt);
    fclose(I);
    // write the SA to file
    if(printca){
        cafile = basename + std::string(".gca");
        if((gCa = fopen(cafile.c_str(), "w")) == nullptr){
            std::cerr << "open() file " + cafile + " failed!" << std::endl;
            exit(-1);
        }
        // write the whole SA in a single operation
        if(fwrite(&CA[0],sizeof(uint_s),n,gCa)!=n){ std::cerr << "gCA file write error... exiting!" << std::endl; exit(-1); }
        // close output file
        fclose(gCa);
    }
    end = std::chrono::steady_clock::now();
    if(arg.verbose){
        std::cout << "Writing the eBWT (and gCA): Elapsed time in seconds: "
        << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
        << " s" << std::endl;
    }
}

// function for constructing and storing the BWT and the conjugate array
// of a text without an end of string character ($)
void compute_bwt_wo_dol(Args arg){ 
    // input text initialization
    std::vector<unsigned char> Text;
    uint_s n=0;
    // initalize output files
    std::string bwtfile, ifile, cafile;
    FILE *Bwt, *I, *Ca;
    // initialize output basename
    std::string basename = arg.outname;
    // initialize boolean values
    bool printca = arg.ca;
    // read input
    load_text(arg.filename.c_str(), Text, n);
    // initialize the conjugate array
    std::vector<uint_s> CA(n,0);
    // initialize the names of the output files
    // compute the conjugate array
    auto start = std::chrono::steady_clock::now();
    cais_bwt(&Text[0], &CA[0], n, alph_size);
    auto end = std::chrono::steady_clock::now();
    if(arg.verbose){
        std::cout << "CA construction: Elapsed time in seconds: "
            << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
            << " s" << std::endl;
    }
    start = std::chrono::steady_clock::now();
    // initialize the names of the output files
    bwtfile = basename + std::string(".bwt");
    ifile = basename + std::string(".i");
    // open output files
    if((Bwt = fopen(bwtfile.c_str(), "w")) == nullptr){
        std::cerr << "open() file " + bwtfile + " failed!" << std::endl;
        exit(-1);
    }
    if((I = fopen(ifile.c_str(), "w")) == nullptr){
        std::cerr << "open() file " + ifile + " failed!" << std::endl;
        exit(-1);
    }
    // write the bwt and the i vector to file
    for(size_t i=0; i<CA.size(); ++i){
        char c;
        // if CA[i] is not the first position of the text
        if( CA[i]>0 ){
            c = Text[CA[i]-1];
        }else{
            c = Text[n-1];
            if( fwrite(&i, sizeof(uint_s), 1, I) !=1 ){ 
                std::cerr << "I file write error... exiting!" << std::endl;
                exit(1); }
        }
        // write a new character in the bwt
        putc(c, Bwt);
    }
    // close output files
    fclose(Bwt);
    fclose(I);
    // write the SA to file
    if(printca){
        cafile = basename + std::string(".ca");
        if((Ca = fopen(cafile.c_str(), "w")) == nullptr){
            std::cerr << "open() file " + cafile + " failed!" << std::endl;
            exit(-1);
        }
        // write the whole SA in a single operation
        if( fwrite(&CA[0], sizeof(uint_s), n, Ca) !=n ){ std::cerr << "CA file write error... exiting!" << std::endl; exit(-1); }
        // close output file
        fclose(Ca);
    }
    end = std::chrono::steady_clock::now();
    if(arg.verbose){
    std::cout << "Write the BWT: Elapsed time in seconds: "
        << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
        << " s" << std::endl;
    }
}

// function for constructing and storing the BBWT and 
// the generalized conjugate array of a text
void compute_bbwt(Args arg){ 
    // input text initialization
    std::vector<unsigned char> Text = {};
    uint_s n=0;
    std::vector<uint_s> onset = {};
    // initialize gap encoded bit-vector 
    sdsl::sd_vector<> b;
    sdsl::sd_vector<>::rank_1_type r_s;
    sdsl::sd_vector<>::select_1_type s_s;
    // initalize output files
    std::string bbwtfile, ifile, cafile;
    FILE *Bbwt, *I, *gCa;
    // initialize output basename
    std::string basename = arg.outname;
    // initialize boolean values
    bool printca = arg.ca;
    // read input
    load_text(arg.filename.c_str(), Text, n);
    std::cout << "size: " << n << std::endl;
    // compute the Lyndon factorization in linear time
    size_t i = 0;
    while (i < n) {
      size_t j = i + 1, k = i;
      while (j < n && Text[k] <= Text[j]) {
          if (Text[k] < Text[j])
              k = i;
          else
              k++;
          j++;
      }
      while (i <= k) {
        onset.push_back(i);
        // j - k cannot be negative
        i += j - k; 
      }
    }
    // add last element in the onset vector
    onset.push_back(n);
    // construct gap-encoded bit vector marking the starting position
    // of each Lyndon factor
    sdsl::sd_vector_builder builder(n+1,onset.size());
    for(auto idx: onset){ builder.set(idx); }
    b = sdsl::sd_vector<>(builder);
    // free onset vector
    onset.clear();
    // initialize the generalized conjugate array
    std::vector<uint_s> CA(n,0);
    // initialize the names of the output files
    // compute the generalized conjugate array
    auto start = std::chrono::steady_clock::now();
    cais(&Text[0], &CA[0], n, alph_size, b);
    auto end = std::chrono::steady_clock::now();
    if(arg.verbose){
    std::cout << "gCA construction: Elapsed time in seconds: "
        << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
        << " s" << std::endl;
    }
    start = std::chrono::steady_clock::now();
    // initialize the names of the output files
    bbwtfile = basename + std::string(".bbwt");
    ifile = basename + std::string(".bbi");
    // open output files
    if((Bbwt = fopen(bbwtfile.c_str(), "w")) == nullptr){
        std::cerr << "open() file " + bbwtfile + " failed!" << std::endl;
        exit(-1);
    }
    if((I = fopen(ifile.c_str(), "w")) == nullptr){
        std::cerr << "open() file " + ifile + " failed!" << std::endl;
        exit(-1);
    }
    // initialize rank and select data structures
    r_s = sdsl::sd_vector<>::rank_1_type(&b);
    s_s = sdsl::sd_vector<>::select_1_type(&b);
    // write the bbwt and the i vector to file
    for(size_t i=0; i<CA.size(); ++i){
        char c;
        // if CA[i] is not the first position of a string
        if(b[CA[i]]!=1){ c = Text[CA[i]-1];  }
        // write the position in the i vector
        else{ 
            c = Text[s_s(r_s(CA[i]+1)+1)-1];
            if(fwrite(&i, sizeof(uint32_t), 1, I) !=1 ){
                std::cerr << "I file write error... exiting!" << std::endl;
                exit(-1);
            }
          }
        // write a new character in the ebwt
        putc(c, Bbwt);
    }
    // close output files
    fclose(Bbwt);
    fclose(I);
    // write the SA to file
    if(printca){
        cafile = basename + std::string(".gca");
        if((gCa = fopen(cafile.c_str(), "w")) == nullptr){
            std::cerr << "open() file " + cafile + " failed!" << std::endl;
            exit(-1);
        }
        // write the whole SA in a single operation
        if( fwrite(&CA[0], sizeof(uint_s), n, gCa) !=n ){ std::cerr << "CA file write error... exiting!" << std::endl; exit(-1); }
        // close output file
        fclose(gCa);
    }
    end = std::chrono::steady_clock::now();
    if(arg.verbose){
        std::cout << "write BBWT: Elapsed time in seconds: "
        << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
        << " s" << std::endl;
    }
}