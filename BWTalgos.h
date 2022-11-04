#ifndef BWTVAR_H  
#define BWTVAR_H

#include <vector>
#include <iostream>
#include <string>
#include <chrono>
#include <unistd.h>

// algo for ca and gca construction 
#include "lib/cais.h"
#include "IOfunc.hpp"

#include "external/malloc_count/malloc_count.h"
#include "external/malloc_count/stack_count.h"

// alphabet size
const int alph_size = 128;

void compute_ebwt(Args arg, bool concat);
void compute_bwt_wo_dol(Args arg);
void compute_bbwt(Args arg);

#endif