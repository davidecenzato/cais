# compilation flags
CXX_FLAGS=-std=c++11 -O3 -Wall -Wextra -pedantic -g
#CXX_FLAGS=-std=c++11 -O3 -Wall -Wextra -pedantic -fsanitize=address -fno-omit-frame-pointer -g
CFLAGS=-O3 -Wall -std=c99 -g
CC=gcc
CCX=g++

# main executables 
EXECS = cais cais64

# targets not producing a file declared phony
.PHONY: all clean

all: $(EXECS)

lib/cais32.o: lib/cais.cpp lib/cais.h
	$(CCX) $(CXX_FLAGS) -c -o $@ $< -lsdsl -ldivsufsort -ldivsufsort64

lib/cais64.o: lib/cais.cpp lib/cais.h
	$(CCX) $(CXX_FLAGS) -c -o $@ $< -lsdsl -ldivsufsort -ldivsufsort64 -DM64

cais: main.cpp IOfunc.hpp BWTalgos.cpp external/malloc_count/malloc_count.o lib/cais32.o 
	$(CXX) $(CXX_FLAGS) -o $@ main.cpp IOfunc.hpp BWTalgos.cpp external/malloc_count/malloc_count.o lib/cais32.o -ldl -lsdsl -ldivsufsort -ldivsufsort64

cais64: main.cpp IOfunc.hpp BWTalgos.cpp external/malloc_count/malloc_count.o lib/cais64.o 
	$(CXX) $(CXX_FLAGS) -o $@ main.cpp IOfunc.hpp BWTalgos.cpp external/malloc_count/malloc_count.o lib/cais64.o -ldl -lsdsl -ldivsufsort -ldivsufsort64 -DM64

clean:
	rm -f $(EXECS) $(EXECS_NT) lib/*.o  