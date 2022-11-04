/*
 * Conjugate Array Induced Sorting (cais) implementation to compute either the generalized conjugate array
 * of a string collection or the conjugate array of a text without dollar.
 * 
 * This code is adapted from https://github.com/kurpicz/saca-bench/blob/master/sa-is/sais.cpp
 * which is the original code of the SA-IS algorithm listed below, and 
 * from https://github.com/felipelouza/gsa-is/blob/master/gsacak.c which is an implementation
 * of the GSACA-K algorithm.
 *
 * void cais(s, CA, n, K, b_s); // compute the generalized conjugate array of a string collection.
 * void cais_bwt(s, CA, n, K);  // compute the conjugate array of a text without dollar.
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
    typedef uint64_t uint_s;
    typedef int64_t int_s;
    #define U_MAX	UINT64_MAX
#else
    typedef uint32_t uint_s;
    typedef int32_t int_s;
    #define U_MAX	UINT32_MAX
#endif
// boolean type
typedef bool bool_t;
// empty value
const uint_s EMPTY = U_MAX;

/**
 * @brief Compute the generalized conjugate array of a string collection using cais algo
 * 
 * @param s input string
 * @param CA generalized conjugate array
 * @param n length of the input string
 * @param K alphabet size
 * @param b_s bitvector of the starting phrases of the parse
 */
void cais(unsigned char *s, uint_s *CA, size_t n, size_t K, sdsl::sd_vector<> &b_s);

/**
 * @brief Compute the conjugate array of a text using cais algo
 * 
 * @param s input string
 * @param CA conjugate array
 * @param n length of the input string
 * @param K alphabet size
 */
void cais_bwt(unsigned char *s, uint_s *CA, size_t n, size_t K);


#endif /* CAIS_H */

