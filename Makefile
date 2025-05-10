# compilation flags
CCX_FLAGS=-std=c++11 -O3 -Wall -Wextra -pedantic -g -I external/sdsl-lite/installed/include -L external/sdsl-lite/installed/lib
CCX=g++

# main executables 
EXECS = cais cais64

# targets not producing a file declared phony
.PHONY: all clean

all: $(EXECS)

lib/cais32.o: lib/cais.cpp lib/cais.h
	$(CCX) $(CCX_FLAGS) -c -o $@ $< -lsdsl -ldivsufsort -ldivsufsort64

lib/cais64.o: lib/cais.cpp lib/cais.h
	$(CCX) $(CCX_FLAGS) -c -o $@ $< -lsdsl -ldivsufsort -ldivsufsort64 -DM64

cais: main.cpp BWTalgos.cpp external/malloc_count/malloc_count.o lib/cais32.o 
	$(CCX) $(CCX_FLAGS) -o $@ main.cpp BWTalgos.cpp external/malloc_count/malloc_count.o lib/cais32.o -ldl -lsdsl -ldivsufsort -ldivsufsort64

cais64: main.cpp BWTalgos.cpp external/malloc_count/malloc_count.o lib/cais64.o 
	$(CCX) $(CCX_FLAGS) -o $@ main.cpp BWTalgos.cpp external/malloc_count/malloc_count.o lib/cais64.o -ldl -lsdsl -ldivsufsort -ldivsufsort64 -DM64

clean:
	rm -f $(EXECS) lib/*.o  