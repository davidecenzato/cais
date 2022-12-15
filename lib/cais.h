/*
 * Conjugate Array Induced Sorting (cais) implementation to compute either the generalized conjugate array
 * of a string collection or the conjugate array of a text without dollar.
 * 
 * This code is adapted from https://github.com/kurpicz/saca-bench/blob/master/sa-is/sais.cpp
 * which is the original code of the SA-IS algorithm listed below, and 
 * from https://github.com/felipelouza/gsa-is/blob/master/gsacak.c which is an implementation
 * of the GSACA-K algorithm.
 *
 * void cais(s, CA, n, K, bv, r_bv, s_bv); // compute the generalized conjugate array of a string collection using a sparse bitvector.
 * void cais_plain(s, CA, n, K, bv, r_bv, s_bv); // compute the generalized conjugate array of a string collection using a plain bitvector.
 * void cais_bvwt(s, CA, n, K);  // compute the conjugate array of a text without dollar.
 *
 */

#ifndef CAIS_H
#define CAIS_H

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <inttypes.h>
#include <string.h>
#include <memory.h>

#include <iostream>
#include <ctime>
#include <vector>
#include <cassert>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>

#ifndef M64
	#define M64 0
#endif

#if M64
    typedef uint64_t uint_t;
    typedef int64_t int_t;
    #define U_MAX	UINT64_MAX
#else
    typedef uint32_t uint_t;
    typedef int32_t int_t;
    #define U_MAX	UINT32_MAX
#endif
// boolean type
typedef bool bool_t;
// empty value
const uint_t EMPTY = U_MAX;

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
	      sdsl::sd_vector<>::rank_1_type& r_bv, sdsl::sd_vector<>::select_1_type& s_bv);

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
	            sdsl::bit_vector::rank_1_type& r_bv, sdsl::bit_vector::select_1_type& s_bv);

/**
 * @brief Compute the conjugate array of a text using cais algo
 * 
 * @param s input string
 * @param CA conjugate array
 * @param n length of the input string
 * @param K alphabet size
 */
void cais_bwt(unsigned char *s, uint_t *CA, size_t n, size_t K);


#endif /* CAIS_H */

