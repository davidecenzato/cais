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
// vector of pairs
typedef std::vector<std::pair<uint_t,uint_t>> pair_vector;
// empty value
const uint_t EMPTY = U_MAX;

/**
 * @brief Compute the generalized conjugate array of a string collection of primitive strings using the cais algorithm
 * 
 * @param s     concatenated input string collection
 * @param CA    generalized conjugate array
 * @param n     concatenation length
 * @param K     alphabet size
 * @param bv    sparse sdsl bitvector storing the string delimiters
 * @param r_bv  rank data structure for sdsl sparse bitvector
 * @param s_bv  select data structure for sdsl sparse bitvector
 * @return      none
 */
void cais_primitive(unsigned char *s, uint_t *CA, size_t n, size_t K, sdsl::sd_vector<> &bv,
                    sdsl::sd_vector<>::rank_1_type& r_bv, sdsl::sd_vector<>::select_1_type& s_bv);

/**
 * @brief Compute the generalized conjugate array of any arbitrary string collection using the cais algorithm
 * 
 * @param s     concatenated input string collection
 * @param CA    generalized conjugate array
 * @param EXP   array containing the string exponents
 * @param n     concatenation length
 * @param K     alphabet size
 * @param bv    sparse sdsl bitvector storing the string delimiters
 * @param r_bv  rank data structure for sdsl sparse bitvector
 * @param s_bv  select data structure for sdsl sparse bitvector
 * @return      none
 */
void cais(unsigned char *s, uint_t *CA, uint_t *EXP, size_t n, size_t K, sdsl::sd_vector<> &bv,
          sdsl::sd_vector<>::rank_1_type& r_bv, sdsl::sd_vector<>::select_1_type& s_bv);

/**
 * @brief Compute the generalized conjugate array of a string collection of primitive strings using the cais algorithm
 * 
 * @param s     concatenated input string collection
 * @param CA    generalized conjugate array
 * @param n     concatenation length
 * @param K     alphabet size
 * @param bv    uncompressed sdsl bitvector storing the string delimiters
 * @param r_bv  rank data structure for sdsl uncompressed bitvector
 * @param s_bv  select data structure for sdsl uncompressed bitvector
 * @return      none
 */
void cais_uncomp_primitive(unsigned char *s, uint_t *CA, size_t n, size_t K, sdsl::bit_vector &bv,
                           sdsl::bit_vector::rank_1_type& r_bv, sdsl::bit_vector::select_1_type& s_bv);

/**
 * @brief Compute the generalized conjugate array of any arbitrary string collection using the cais algorithm
 * 
 * @param s     concatenated input string collection
 * @param CA    generalized conjugate array
 * @param EXP   array containing the string exponents
 * @param n     concatenation length
 * @param K     alphabet size
 * @param bv    uncompressed sdsl bitvector storing the string delimiters
 * @param r_bv  rank data structure for sdsl uncompressed bitvector
 * @param s_bv  select data structure for sdsl uncompressed bitvector
 * @return      None
 */
void cais_uncomp(unsigned char *s, uint_t *CA, uint_t *EXP, size_t n, size_t K, sdsl::bit_vector &bv,
                 sdsl::bit_vector::rank_1_type& r_bv, sdsl::bit_vector::select_1_type& s_bv);

/**
 * @brief Compute the conjugate array of one string using the cais algorithm
 * 
 * @param s  input string
 * @param CA conjugate array
 * @param n  length of the input string
 * @param K  alphabet size
 * @return   none
 */
void cais_bwt(unsigned char *s, uint_t *CA, size_t n, size_t K);


#endif /* CAIS_H */