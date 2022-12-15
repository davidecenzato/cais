/*
 * Conjugate Array Induced Sorting (cais) implementation to compute either the generalized conjugate array
 * of a string collection or the conjugate array of a text without dollar.
 * 
 * This code is adapted from https://github.com/kurpicz/saca-bench/blob/master/sa-is/sais.cpp
 * which is the original code of the SA-IS algorithm listed below, and 
 * from https://github.com/felipelouza/gsa-is/blob/master/gsacak.c which is an implementation
 * of the GSACA-K algorithm.
 *
 */

#include "cais.h"

// get the character function
#define chr(i) (cs==sizeof(uint_t)?((uint_t *)s)[i]:((unsigned char *)s)[i])

// compute the head or end of each bucket
void getBuckets(uint_t *s, uint_t *bkt, size_t n, size_t K, size_t cs, bool end) { 
  size_t i, sum=0;
  for(i=0; i<K; i++) bkt[i]=0; // clear all buckets
  for(i=0; i<n; i++) bkt[chr(i)]++;  // compute the size of each bucket
  for(i=0; i<K; i++) { sum+=bkt[i]; bkt[i]= end ? sum-1 : sum-bkt[i]; }
}

// induce L suffixes for cais algorithm using sdsl sparse bitvector
void induceL(uint_t *CA, uint_t *s, uint_t *bkt, sdsl::sd_vector<>::rank_1_type& r_bv, 
             sdsl::sd_vector<>::select_1_type& s_bv, sdsl::sd_vector<>& bv, size_t n, size_t K, 
             size_t cs, bool phase) { 
    
    size_t i, j, m, rank;
    getBuckets(s, bkt, n, K, cs, false); // find heads of buckets

    for(i=0; i<n; ++i){
        if(CA[i]!=EMPTY) {
            j=CA[i];
            if(bv[j]==1){
                rank=r_bv(j+1);m=s_bv(rank+1)-1;
                if(chr(m) >= chr(j)){ CA[bkt[chr(m)]++]=m; if(phase){ CA[i] = EMPTY; } }
            }else{
                if(chr(j-1) >= chr(j)){ CA[bkt[chr(j-1)]++]=j-1; if(phase){ CA[i] = EMPTY; } } 
            }
        }
    }
}

// induce L suffixes for cais algorithm using sdsl plain bitvector
void induceL_plain(uint_t *CA, uint_t *s, uint_t *bkt, sdsl::bit_vector::rank_1_type& r_bv, 
             sdsl::bit_vector::select_1_type& s_bv, sdsl::bit_vector& bv, size_t n, size_t K, 
             size_t cs, bool phase) { 
    
    size_t i, j, m, rank;
    getBuckets(s, bkt, n, K, cs, false); // find heads of buckets

    for(i=0; i<n; ++i){
        if(CA[i]!=EMPTY) {
            j=CA[i];
            if(bv[j]==1){
                rank=r_bv(j+1);m=s_bv(rank+1)-1;
                if(chr(m) >= chr(j)){ CA[bkt[chr(m)]++]=m; if(phase){ CA[i] = EMPTY; } }
            }else{
                if(chr(j-1) >= chr(j)){ CA[bkt[chr(j-1)]++]=j-1; if(phase){ CA[i] = EMPTY; } } 
            }
        }
    }
}

// induce L suffixes for cais bwt algorithm
void induceLbwt(uint_t *CA, uint_t *s, uint_t *bkt, size_t n, size_t K, 
                size_t cs, bool phase) { 
    
    size_t i, j;
    getBuckets(s, bkt, n, K, cs, false); // find heads of buckets

    for(i=0; i<n; i++){
        if(CA[i]!=EMPTY) {
            j=CA[i];
            if(j==0){
                if(chr(n-1) >= chr(0)){ CA[bkt[chr(n-1)]++]=n-1; if(phase){ CA[i] = EMPTY; } }
            }else{
                if(chr(j-1) >= chr(j)){ CA[bkt[chr(j-1)]++]=j-1; if(phase){ CA[i] = EMPTY; } } 
            }
        }
    }
}

// induce S suffixes for cais algorithm using sdsl sparse bitvector
void induceS(uint_t *CA, uint_t *s, uint_t *bkt, std::vector<uint_t>& singletons, sdsl::sd_vector<>::rank_1_type& r_bv,
             sdsl::sd_vector<>::select_1_type& s_bv, sdsl::sd_vector<>& bv, size_t n, size_t K, size_t cs, bool phase) { 
    
    size_t i, j, m, rank;
    getBuckets(s, bkt, n, K, cs, true); // find ends of buckets

    for(i=0; i<n; ++i){
        size_t ni = n-i-1;
        if(CA[ni]!=EMPTY) {
            j=CA[ni];
            if(bv[j]==1){
                rank=r_bv(j+1);m=s_bv(rank+1)-1;
                if(chr(m) <= chr(j) && bkt[chr(m)]<ni){ CA[bkt[chr(m)]--]=m; if(phase){ CA[ni] = EMPTY; } }    
             }else{
                if(chr(j-1) <= chr(j) && bkt[chr(j-1)]<ni){ CA[bkt[chr(j-1)]--]=j-1; if(phase){ CA[ni] = EMPTY; } } 
            }  
        }   
    }
  // insert singletons
  if(!phase){ 
    // scan all singletons
    size_t sz = singletons.size();
    for(i = 0;i < sz; ++i){
      //uint_t curr = singletons[sz-i-1];
      uint_t curr = singletons[i];
      size_t slen = curr - s_bv(r_bv(curr+1)) + 1;
      for(j = 0; j < slen; ++j){ CA[bkt[chr(curr-j)]--]=curr-j; } 
    }
  }
}

// induce S suffixes for cais algorithm using sdsl plain bitvector
void induceS_plain(uint_t *CA, uint_t *s, uint_t *bkt, std::vector<uint_t>& singletons, sdsl::bit_vector::rank_1_type& r_bv,
             sdsl::bit_vector::select_1_type& s_bv, sdsl::bit_vector& bv, size_t n, size_t K, size_t cs, bool phase) { 
    
    size_t i, j, m, rank;
    getBuckets(s, bkt, n, K, cs, true); // find ends of buckets

    for(i=0; i<n; ++i){
        size_t ni = n-i-1;
        if(CA[ni]!=EMPTY) {
            j=CA[ni];
            if(bv[j]==1){
                rank=r_bv(j+1);m=s_bv(rank+1)-1;
                if(chr(m) <= chr(j) && bkt[chr(m)]<ni){ CA[bkt[chr(m)]--]=m; if(phase){ CA[ni] = EMPTY; } }    
             }else{
                if(chr(j-1) <= chr(j) && bkt[chr(j-1)]<ni){ CA[bkt[chr(j-1)]--]=j-1; if(phase){ CA[ni] = EMPTY; } } 
            }  
        }   
    }
  // insert singletons
  if(!phase){ 
    // scan all singletons
    size_t sz = singletons.size();
    for(i = 0;i < sz; ++i){
      //uint_t curr = singletons[sz-i-1];
      uint_t curr = singletons[i];
      size_t slen = curr - s_bv(r_bv(curr+1)) + 1;
      for(j = 0; j < slen; ++j){ CA[bkt[chr(curr-j)]--]=curr-j; } 
    }
  }
}

// induce S suffixes for cais bwt algorithm
void induceSbwt(uint_t *CA, uint_t *s, uint_t *bkt, size_t n, size_t K, size_t cs, bool phase) { 
    
    size_t i, j;
    getBuckets(s, bkt, n, K, cs, true); // find ends of buckets
    for(i=0; i<n; ++i){
        size_t ni = n-i-1;
        if(CA[ni]!=EMPTY) {
            j=CA[ni];
            if(j==0){
                if(chr(n-1) <= chr(0) && bkt[chr(n-1)]<ni){ CA[bkt[chr(n-1)]--]=n-1; if(phase){ CA[ni] = EMPTY; } }    
             }else{
                if(chr(j-1) <= chr(j) && bkt[chr(j-1)]<ni){ CA[bkt[chr(j-1)]--]=j-1; if(phase){ CA[ni] = EMPTY; } } 
            }  
        }   
    }
}

// compute length of a LMS substring for cais algorithm
uint_t LMSlength(uint_t *s, uint_t sb, uint_t eb, uint_t x, size_t cs) {
  // return the length of the LMS substring
  uint_t len=1, i=1;  
  uint_t prev = x, pos = 0;
  if(x == eb){ pos = sb; }else{ pos = x+1; }
  // S suffixes
  while(true) {
    if(chr(pos)<chr(prev)){ break; }
    ++i;
    if(pos == eb) { pos = sb; } else{ ++pos;  }
    if(prev == eb){ prev = sb; }else{ ++prev; }
  }  
  // L suffixes
  while(true) {
    if(chr(pos)>chr(prev)){ break; }
    if(chr(pos)<chr(prev)){ len=i; }
    ++i;
    if(pos == eb) { pos = sb; } else{ ++pos;  }
    if(prev == eb){ prev = sb; }else{ ++prev; }
  }
  // return the length
  return len+1;
}

/** @brief compute the generalized conjugate array of string s[0..n-1] using CAIS algo
 *
 *  @param s     input string 
 *  @param CA    generalized conjugate array 
 *  @param n     string length
 *  @param K     alphabet size
 *  @param cs    integer size
 *  @param bv   sparse bitvector for string delimiters
 *  @param r_bv   support for rank queries
 *  @param s_bv   support for select queries
 *  @return      None
 */
void CAIS(uint_t *s, uint_t *CA, size_t n, size_t K, size_t cs, sdsl::sd_vector<> &bv,
          sdsl::sd_vector<>::rank_1_type& r_bv, sdsl::sd_vector<>::select_1_type& s_bv) {
  // initialize variables
  size_t i, j, m, nseq;
  size_t sb, eb, fl, len; 
  // singletons vector
  std::vector<uint_t> sn;
  // no. sequences
  nseq = r_bv(n);
  // stage 1: reduce the problem by at least 1/2
  uint_t *bkt = (uint_t *)malloc(sizeof(uint_t)*K); // bucket counters
  
  for(i=0; i<n; i++) CA[i]=EMPTY; // initialize CA values to -1
  getBuckets(s, bkt, n, K, cs, true); // find ends of buckets
  
  //initialize onset vector
  std::vector<uint_t> onset; size_t st = 0;
  onset.reserve(nseq+1); onset.push_back(st);
  // place S* suffixes in their buckets
  for(i=0; i<nseq; ++i){
      sb = s_bv(i+1), eb = s_bv(i+2)-1, fl = eb+1, len=eb-sb+1;
      assert(len > 0);
      bool type;
      if(len > 1){
        // find first type
        for(j=eb; j>sb; --j){
            if(chr(j)!=chr(j-1)){
              if(chr(j)>chr(j-1)){fl=j; type=0; break;}
              else{fl=j-1; type=0; break;}
            } 
        }
        // set other types
        if(fl<eb+1){
          // Classify the S* suffixes
          uint_t pos, prev; prev=fl; 
          if(prev==sb){ pos = eb; } else{ pos = fl-1; }
          while(pos != fl){
              if(chr(pos) > chr(prev)){ if(type){ CA[bkt[chr(prev)]--] = prev; ++st; type = 0; } }
              else{ 
                  if(chr(pos)<chr(prev)){ type = 1; }
              }
              // Skip to next chars
              if(prev==sb){ prev = eb; } else{ --prev; }
              if(pos==sb) { pos = eb;  } else{ --pos;  }
           }
          // check the last suffix
          if(chr(pos) > chr(prev)){ if(type){ CA[bkt[chr(prev)]--] = prev; ++st; type = 0; } }
          onset.push_back(st);
        }
        //onset.push_back(st); 
     }
  } 
  
  // Induce L ans S suffixes to sort the circular LMS substrings
  induceL(CA, s, bkt, r_bv, s_bv, bv, n, K, cs, true);
  induceS(CA, s, bkt, sn, r_bv, s_bv, bv, n, K, cs, true); 
  
  free(bkt); // free bucket vector
  
  // compact all the sorted substrings into the first n1 items of s
  size_t n1=0;
  for(i=0; i<n; ++i){
      uint_t pos=CA[i];
      if(pos != EMPTY){ CA[n1++]=pos; }
  }
  
  // Init the name array buffer
  for(i=n1; i<n; ++i) CA[i]=EMPTY;
  
  // find the lexicographic names of all LMS substrings
  uint_t name=1;
  if(n1 > 0){
    // initialize names array
    uint_t *names = (uint_t *)malloc(sizeof(uint_t)*n); for(i=0; i<n; i++) names[i]=EMPTY;
    // insert the first LMS substring
    uint_t pos = CA[0];
    names[pos]=name-1;
    size_t rank = r_bv(pos+1), sb = s_bv(rank), eb = s_bv(rank+1)-1;
    uint_t pre_len = LMSlength(s, sb, eb, pos, cs); 
    uint_t prev = pos;  size_t pre_sb = sb, pre_eb = eb;
    // for all S* suffixes
    for(i=1; i<n1; ++i) {
          pos=CA[i]; bool diff=false;
          rank = r_bv(pos+1); sb = s_bv(rank); eb = s_bv(rank+1)-1;
          uint_t len = LMSlength(s, sb, eb, pos, cs); 
          // if the LMS length are different skip and increase name counter
          if(len != pre_len){ diff = true; }
          else{
              // if same length compare all the characters till you find a mismatch
              uint_t tp = pos, pre_tp = prev;
              for(j=0;j<len;++j){ 
                  if(chr(pre_tp)!=chr(tp)){ diff=true; break; }
                  if(tp==eb){ tp=sb; }else{ ++tp;}
                  if(pre_tp==pre_eb){ pre_tp=pre_sb;}else{ ++pre_tp;}
              }
          }
        // name the LMS substring
        if(diff){ ++name; prev=pos; pre_len = len; pre_sb = sb; pre_eb = eb; }
        names[pos]=name-1;
    }
    // compact names array
    j = n-1;
    for(i=0; i<n; ++i){
        if(names[n-i-1]!=EMPTY){CA[j--]=names[n-i-1];}
    }
    // free names array
    free(names);
  }
  
  // s1 is done now
  uint_t *CA1=CA, *s1=CA+n-n1;  

  // stage 2: solve the reduced problem
  // recurse if names are not unique yet
  if(name<n1) {
      // compute the new bit vector
      sdsl::sd_vector_builder builder(n1+1,onset.size());
      for(auto idx: onset){ builder.set(idx); }
      sdsl::sd_vector<> nbv = sdsl::sd_vector<>(builder);
      // free onset vector
      onset.clear(); 
      // compute support for rank and select queries 
      sdsl::sd_vector<>::rank_1_type nr_bv = sdsl::sd_vector<>::rank_1_type(&nbv);
      sdsl::sd_vector<>::select_1_type ns_bv = sdsl::sd_vector<>::select_1_type(&nbv);
      // recursive call
      CAIS((uint_t*)s1, CA1, n1, name, sizeof(uint_t), nbv, nr_bv, ns_bv);
  } else { // stop the recursion, generate the suffix array of s1 directly
      for(i=0; i<n1; i++){ CA1[s1[i]] = i; }
        // free last onset vector
        onset.clear();
  }
  // stage 3: induce the result for the original problem
  bkt = (uint_t *)malloc(sizeof(uint_t)*K); // bucket counters
  
  // put all left-most S characters into their buckets
  j=n1-1;
  for(i=0; i<nseq; ++i){
      size_t ni = nseq-i; bool type, mismatch;
      sb=s_bv(ni); eb=s_bv(ni+1)-1; len=eb-sb+1;
      if(len==1){ sn.push_back(sb); } // fill singletons vector 
      else{
          // if it is not a singleton find type of the first int
          type = 0, mismatch = 0;
          if(chr(sb)!=chr(sb+1)){ if(chr(sb)<chr(sb+1)){ type = 1; } mismatch = 1; }
          else{
              m = 1;
              while(m < len-1){
                  if(chr(sb+m)!=chr(sb+1+m)){ type=(chr(sb+m)>chr(sb+1+m))?0:1; mismatch = 1; break; }
                  ++m;
              }
          }
          // type of last int
          if(chr(eb)!=chr(sb)){ type = (chr(eb)>chr(sb))?0:1; }
          // if the sequence is a power
          if(!mismatch){ sn.push_back(eb); }
          else{
            // find all S* suffixes
            for(m=0; m<len-1; ++m){
                size_t nx = eb-m-1;
                if(chr(nx)>chr(nx+1)){ if(type){ s1[j--]=nx+1; type=0; } }
                else{ if(chr(nx)<chr(nx+1)){ type=1; } }
            }
            // check if first suffix is S*
            if(chr(sb)<chr(eb)){ if(type){ s1[j--]=sb; } }
          }
      }
  }
  
  // we must have the same number of S* of stage 1
  assert(j+1==0);
    
  for(i=0; i<n1; ++i) {CA1[i]=s1[CA1[i]];} // get index in s1
   
  if(n1>0){ // if at least one s* is not a singleton
      for(i=n1; i<n; ++i) CA[i]=EMPTY; // init CA[n1..n-1]
      getBuckets(s, bkt, n, K, cs, true); // find ends of buckets
  }
  
  // insert S* suffixes at the end of the buckets
  for(i=0; i < n1; ++i){
      j=CA[n1-i-1]; CA[n1-i-1]=EMPTY;
      CA[bkt[chr(j)]--]=j;
  }
  
  // induce the L and S conjugates to compute the final CA of each level
  induceL(CA, s, bkt, r_bv, s_bv, bv, n, K, cs, false);
  induceS(CA, s, bkt, sn, r_bv, s_bv, bv, n, K, cs, false); 
  
  // free bucket vector
  free(bkt); 
}

/** @brief compute the generalized conjugate array of string s[0..n-1] using CAIS algo
 *
 *  @param s     input string 
 *  @param CA    generalized conjugate array 
 *  @param n     string length
 *  @param K     alphabet size
 *  @param cs    integer size
 *  @param bv   bitvector for string delimiters
 *  @param r_bv   support for rank queries
 *  @param s_bv   support for select queries
 *  @return      None
 */
void CAIS_P(uint_t *s, uint_t *CA, size_t n, size_t K, size_t cs, sdsl::bit_vector &bv,
            sdsl::bit_vector::rank_1_type& r_bv, sdsl::bit_vector::select_1_type& s_bv) {
  // initialize variables
  size_t i, j, m, nseq;
  size_t sb, eb, fl, len; 
  // singletons vector
  std::vector<uint_t> sn;
  // no. sequences
  nseq = r_bv(n);
  // stage 1: reduce the problem by at least 1/2
  uint_t *bkt = (uint_t *)malloc(sizeof(uint_t)*K); // bucket counters
  
  for(i=0; i<n; i++) CA[i]=EMPTY; // initialize CA values to -1
  getBuckets(s, bkt, n, K, cs, true); // find ends of buckets
  
  //initialize onset vector
  std::vector<uint_t> onset; size_t st = 0;
  onset.reserve(nseq+1); onset.push_back(st);
  // place S* suffixes in their buckets
  for(i=0; i<nseq; ++i){
      sb = s_bv(i+1), eb = s_bv(i+2)-1, fl = eb+1, len=eb-sb+1;
      assert(len > 0);
      bool type;
      if(len > 1){
        // find first type
        for(j=eb; j>sb; --j){
            if(chr(j)!=chr(j-1)){
              if(chr(j)>chr(j-1)){fl=j; type=0; break;}
              else{fl=j-1; type=0; break;}
            } 
        }
        // set other types
        if(fl<eb+1){
          // Classify the S* suffixes
          uint_t pos, prev; prev=fl; 
          if(prev==sb){ pos = eb; } else{ pos = fl-1; }
          while(pos != fl){
              if(chr(pos) > chr(prev)){ if(type){ CA[bkt[chr(prev)]--] = prev; ++st; type = 0; } }
              else{ 
                  if(chr(pos)<chr(prev)){ type = 1; }
              }
              // Skip to next chars
              if(prev==sb){ prev = eb; } else{ --prev; }
              if(pos==sb) { pos = eb;  } else{ --pos;  }
           }
          // check the last suffix
          if(chr(pos) > chr(prev)){ if(type){ CA[bkt[chr(prev)]--] = prev; ++st; type = 0; } }
          onset.push_back(st);
        }
        //onset.push_back(st); 
     }
  } 
  
  // Induce L ans S suffixes to sort the circular LMS substrings
  induceL_plain(CA, s, bkt, r_bv, s_bv, bv, n, K, cs, true);
  induceS_plain(CA, s, bkt, sn, r_bv, s_bv, bv, n, K, cs, true); 
  
  free(bkt); // free bucket vector
  
  // compact all the sorted substrings into the first n1 items of s
  size_t n1=0;
  for(i=0; i<n; ++i){
      uint_t pos=CA[i];
      if(pos != EMPTY){ CA[n1++]=pos; }
  }
  
  // Init the name array buffer
  for(i=n1; i<n; ++i) CA[i]=EMPTY;
  
  // find the lexicographic names of all LMS substrings
  uint_t name=1;
  if(n1 > 0){
    // initialize names array
    uint_t *names = (uint_t *)malloc(sizeof(uint_t)*n); for(i=0; i<n; i++) names[i]=EMPTY;
    // insert the first LMS substring
    uint_t pos = CA[0];
    names[pos]=name-1;
    size_t rank = r_bv(pos+1), sb = s_bv(rank), eb = s_bv(rank+1)-1;
    uint_t pre_len = LMSlength(s, sb, eb, pos, cs); 
    uint_t prev = pos;  size_t pre_sb = sb, pre_eb = eb;
    // for all S* suffixes
    for(i=1; i<n1; ++i) {
          pos=CA[i]; bool diff=false;
          rank = r_bv(pos+1); sb = s_bv(rank); eb = s_bv(rank+1)-1;
          uint_t len = LMSlength(s, sb, eb, pos, cs); 
          // if the LMS length are different skip and increase name counter
          if(len != pre_len){ diff = true; }
          else{
              // if same length compare all the characters till you find a mismatch
              uint_t tp = pos, pre_tp = prev;
              for(j=0;j<len;++j){ 
                  if(chr(pre_tp)!=chr(tp)){ diff=true; break; }
                  if(tp==eb){ tp=sb; }else{ ++tp;}
                  if(pre_tp==pre_eb){ pre_tp=pre_sb;}else{ ++pre_tp;}
              }
          }
        // name the LMS substring
        if(diff){ ++name; prev=pos; pre_len = len; pre_sb = sb; pre_eb = eb; }
        names[pos]=name-1;
    }
    // compact names array
    j = n-1;
    for(i=0; i<n; ++i){
        if(names[n-i-1]!=EMPTY){CA[j--]=names[n-i-1];}
    }
    // free names array
    free(names);
  }
  
  // s1 is done now
  uint_t *CA1=CA, *s1=CA+n-n1;  

  // stage 2: solve the reduced problem
  // recurse if names are not unique yet
  if(name<n1) {
      // construct the new bit vector
      sdsl::bit_vector nbv(n1+1,0);
      for(auto idx: onset){ nbv[idx] = 1; }
      // free onset
      onset.clear();
      // compute support for rank and select queries 
      sdsl::bit_vector::rank_1_type nr_bv = sdsl::bit_vector::rank_1_type(&nbv);
      sdsl::bit_vector::select_1_type ns_bv = sdsl::bit_vector::select_1_type(&nbv);
      // recursive call
      CAIS_P((uint_t*)s1, CA1, n1, name, sizeof(uint_t), nbv, nr_bv, ns_bv);
  } else { // stop the recursion, generate the suffix array of s1 directly
      for(i=0; i<n1; i++){ CA1[s1[i]] = i; }
      // free last onset vector
      onset.clear();
  }
  // stage 3: induce the result for the original problem
  bkt = (uint_t *)malloc(sizeof(uint_t)*K); // bucket counters
  
  // put all left-most S characters into their buckets
  j=n1-1;
  for(i=0; i<nseq; ++i){
      size_t ni = nseq-i; bool type, mismatch;
      sb=s_bv(ni); eb=s_bv(ni+1)-1; len=eb-sb+1;
      if(len==1){ sn.push_back(sb); } // fill singletons vector 
      else{
          // if it is not a singleton find type of the first int
          type = 0, mismatch = 0;
          if(chr(sb)!=chr(sb+1)){ if(chr(sb)<chr(sb+1)){ type = 1; } mismatch = 1; }
          else{
              m = 1;
              while(m < len-1){
                  if(chr(sb+m)!=chr(sb+1+m)){ type=(chr(sb+m)>chr(sb+1+m))?0:1; mismatch = 1; break; }
                  ++m;
              }
          }
          // type of last int
          if(chr(eb)!=chr(sb)){ type = (chr(eb)>chr(sb))?0:1; }
          // if the sequence is a power
          if(!mismatch){ sn.push_back(eb); }
          else{
            // find all S* suffixes
            for(m=0; m<len-1; ++m){
                size_t nx = eb-m-1;
                if(chr(nx)>chr(nx+1)){ if(type){ s1[j--]=nx+1; type=0; } }
                else{ if(chr(nx)<chr(nx+1)){ type=1; } }
            }
            // check if first suffix is S*
            if(chr(sb)<chr(eb)){ if(type){ s1[j--]=sb; } }
          }
      }
  }
  
  // we must have the same number of S* of stage 1
  assert(j+1==0);
    
  for(i=0; i<n1; ++i) {CA1[i]=s1[CA1[i]];} // get index in s1
   
  if(n1>0){ // if at least one s* is not a singleton
      for(i=n1; i<n; ++i) CA[i]=EMPTY; // init CA[n1..n-1]
      getBuckets(s, bkt, n, K, cs, true); // find ends of buckets
  }
  
  // insert S* suffixes at the end of the buckets
  for(i=0; i < n1; ++i){
      j=CA[n1-i-1]; CA[n1-i-1]=EMPTY;
      CA[bkt[chr(j)]--]=j;
  }
  
  // induce the L and S conjugates to compute the final CA of each level
  induceL_plain(CA, s, bkt, r_bv, s_bv, bv, n, K, cs, false);
  induceS_plain(CA, s, bkt, sn, r_bv, s_bv, bv, n, K, cs, false); 
  
  // free bucket vector
  free(bkt); 
}

/** @brief compute the conjugate array of string s[0..n-1] using CAIS bwt algo
 *
 *  @param s     input string 
 *  @param CA    conjugate array 
 *  @param n     string length
 *  @param K     alphabet size
 *  @param cs    integer size
 *  @return      None
 */
void CAISbwt(uint_t *s, uint_t *CA, size_t n, size_t K, size_t cs) {
    size_t i, j, m, fl;
    // stage 1: reduce the problem by at least 1/2
    uint_t *bkt = (uint_t *)malloc(sizeof(uint_t)*K); // bucket counters
    for(i=0; i<n; ++i) CA[i]=EMPTY; // initialize CA values to -1
    getBuckets(s, bkt, n, K, cs, true); // find ends of buckets
    bool type; fl = n; // 0 for L and 1 for S
    for(i=0;i<n-1;++i){
        j = n-i-1;
        if(chr(j) != chr(j-1)){
            if(chr(j) > chr(j-1)){fl=j;}
            else{fl=j-1;}
            type=0; break;
        }
    }
    if(fl != n){ // if the string is not a power
        // search types
        for(i=0;i<fl;++i){
            j = fl-i;
            if(chr(j-1) > chr(j)){
                if(type){ CA[bkt[chr(j)]--] = j; type = 0; } 
            }
            else{ if(chr(j-1)<chr(j)){ type = 1; } }
        }
        // set type of first character of the string
        if(chr(n-1) > chr(0)){ if(type){ CA[bkt[chr(0)]--] = 0; type = 0; } }
        else{ if(chr(n-1)<chr(0)){ type = 1; } }
        // scan all other characters
        if(fl!=n-1){  
            j = n-1;
            while(j != fl){
                if(chr(j-1) > chr(j)){
                    if(type){ CA[bkt[chr(j)]--] = j; type = 0; } 
                    }
                else{ if(chr(j-1)<chr(j)){ type = 1; } }
                j--;
            }
          }
    
        // Induce L ans S suffixes to sort the circular LMS substrings
        induceLbwt(CA, s, bkt, n, K, cs, true);
        induceSbwt(CA, s, bkt, n, K, cs, true);
        
        free(bkt); // free bucket vector
        // compact all the sorted substrings into the first n1 items of s
        size_t n1=0;
        for(i=0; i<n; ++i){
            uint_t pos=CA[i];
            if(pos != EMPTY){ CA[n1++]=pos; }
        }
        
        // Init the name array buffer
        for(i=n1; i<n; ++i) CA[i]=EMPTY;
        
        // find the lexicographic names of all LMS substrings
        uint_t name=1;
        if(n1 > 0){
            // insert the first LMS substring
            uint_t pos = CA[0];
            CA[n1+((pos%2==0)?(pos/2):(pos-1)/2)] = name-1;
            uint_t pre_len = LMSlength(s, 0, n-1, pos, cs); 
            uint_t prev = pos;  
            // for all S* suffixes
            for(i=1; i<n1; i++) {
                pos=CA[i]; bool diff=false;
                uint_t len = LMSlength(s, 0, n-1, pos, cs); 
                // if the LMS length are different skip and increase name counter
                if(len != pre_len){ diff = true; }
                else{
                    // if same length compare all the characters till you find a mismatch
                    uint_t tp = pos, pre_tp = prev;
                    for(j=0;j<len;++j){ 
                        if(chr(pre_tp)!=chr(tp)){ diff=true; break; }
                        if(tp==n-1){ tp=0; }else{ ++tp; }
                        if(pre_tp==n-1){ pre_tp=0; }else{ ++pre_tp; }
                    }
                }   
                // name the LMS substring
                if(diff){ ++name; prev=pos; pre_len = len; }
                CA[n1+((pos%2==0)?(pos/2):(pos-1)/2)] = name-1;
            }
        // compact names array
        j = n-1;
        for(i=n-1; i>=n1; --i){
            if(CA[i]!=EMPTY){CA[j--]=CA[i];}
            }
        }
  
        // s1 is done now
        uint_t *CA1=CA, *s1=CA+n-n1;
        
        // stage 2: solve the reduced problem
        // recurse if names are not unique yet
        if(name<n1) {
            CAISbwt((uint_t*)s1, CA1, n1, name, sizeof(uint_t));
        } else { // stop the recursion, generate the suffix array of s1 directly
            for(i=0; i<n1; ++i){ CA1[s1[i]] = i; }
        }
       
        // stage 3: induce the result for the original problem
        bkt = (uint_t *)malloc(sizeof(uint_t)*K); // bucket counters
        
        // put all left-most S characters into their buckets
        j=n1-1; type = 0;
        // find type of char at index 0
        if(chr(0)!=chr(1)){ if(chr(0)<chr(1)){ type = 1; } }
        else{
            m = 1;
            while(m < n-1){
                if(chr(m)!=chr(m+1)){ type=(chr(m)>chr(m+1))?0:1; break; }
                ++m;
            }
        }
        // find type of last char
        if(chr(n-1)!=chr(0)){ type = (chr(n-1)>chr(0))?0:1; }
        // find all S* suffixes
        for(m=0; m<n-1; ++m){
            size_t nx = n-m-2; //n-1-m-1
            if(chr(nx)>chr(nx+1)){ if(type){ s1[j--]=nx+1; type=0; } }
            else{ 
                if(chr(nx)<chr(nx+1)){ type=1; }
            }
        }
        // check if first suffix is S*
        if(chr(0)<chr(n-1)){ if(type){ s1[j--]=0; } }
        
        // we must have the same number of S* of stage 1
        assert(j+1==0);
        
        for(i=0; i<n1; ++i) {CA1[i]=s1[CA1[i]];} // get index in s1
      
        for(i=n1; i<n; ++i) CA[i]=EMPTY; // init CA[n1..n-1]
        getBuckets(s, bkt, n, K, cs, true); // find ends of buckets
        
        // insert S* suffixes at the end of the buckets
        for(i=0; i < n1; ++i){
            j=CA[n1-i-1]; CA[n1-i-1]=EMPTY;
            CA[bkt[chr(j)]--]=j;
        }
        
        // induce the L and S suffixes to compute the final CA of each level
        induceLbwt(CA, s, bkt, n, K, cs, false);
        induceSbwt(CA, s, bkt, n, K, cs, false); 

        // bucket vector
        free(bkt); 

    }else{
        // if the word is a power
        uint_t *CA1=CA;
        for(i=0; i<n; ++i){ CA1[i] = i; }
        // bucket vector
        free(bkt);
    }
}

/**
 * @brief Compute the generalized conjugate array of a string collection using cais algo
 * 
 * @param s input string
 * @param CA generalized conjugate array
 * @param n length of the input string
 * @param K alphabet size
 * @param bv sparse sdsl bitvector of the starting phrases of the parse
 * @param r_bv rank data structure for sdsl sparse bitvector
 * @param s_bv select data structure for sdsl sparse bitvector
 */
void cais(unsigned char *s, uint_t *CA, size_t n, size_t K, sdsl::sd_vector<> &bv,
          sdsl::sd_vector<>::rank_1_type& r_bv, sdsl::sd_vector<>::select_1_type& s_bv){
    if((s == nullptr) || (CA == nullptr) || (n == 0)) {std::cerr << "Empty input given." << std::endl; exit(1);}
    CAIS((uint_t *)s, (uint_t *)CA, n, K, sizeof(unsigned char), bv, r_bv, s_bv);
}

/**
 * @brief Compute the generalized conjugate array of a string collection using cais algo
 * 
 * @param s input string
 * @param CA generalized conjugate array
 * @param n length of the input string
 * @param K alphabet size
 * @param bv plain sdsl bitvector of the starting phrases of the parse
 * @param r_bv rank data structure for sdsl plain bitvector
 * @param s_bv select data structure for sdsl plain bitvector
 */
void cais_plain(unsigned char *s, uint_t *CA, size_t n, size_t K, sdsl::bit_vector &bv,
                sdsl::bit_vector::rank_1_type& r_bv, sdsl::bit_vector::select_1_type& s_bv){
    if((s == nullptr) || (CA == nullptr) || (n == 0)) {std::cerr << "Empty input given." << std::endl; exit(1);}
    CAIS_P((uint_t *)s, (uint_t *)CA, n, K, sizeof(unsigned char), bv, r_bv, s_bv);
}

/**
 * @brief Compute the conjugate array of a text using cais algo
 * 
 * @param s input string
 * @param CA conjugate array
 * @param n length of the input string
 * @param K alphabet size
 * @param bv bitvector for phrase delimiters
 */
void cais_bwt(unsigned char *s, uint_t *CA, size_t n, size_t K){
    if((s == nullptr) || (CA == nullptr) || (n == 0)) {std::cerr << "Empty input given." << std::endl; exit(-1);}
    CAISbwt((uint_t *)s, (uint_t *)CA, n, K, sizeof(unsigned char));
}

