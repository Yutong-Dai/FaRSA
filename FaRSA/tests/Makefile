# Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
#
# This code is published under the ??? License.
#
# Author(s) : Frank E. Curtis, Daniel P. Robinson

# C++ compiler
CXX = g++

# C++ compiler flags
CXXFLAGS = -g -Wall -std=c++11

# Set sources, etc.
headers = $(wildcard *.hpp)
sources = $(wildcard *.cpp)
objects = $(sources:.cpp=.o)
depends = $(sources:.cpp=.d)

# Executable(s)
EXES = $(sources:.cpp=)

# Libraries
FaRSALIB = $(FARSADIR)/FaRSA/src/libFaRSA.a

# Includes
INCLUDES = -I $(FARSADIR)/FaRSA/src

# Rule for all
all: $(EXES)

# Rule for executable
$(EXES): % : %.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(FaRSALIB) -L $(LAPACKDIR) -ldl -lblas -llapack

# Dependencies for executable
$(EXES): $(FARSALIB)

# Rules for libraries
$(FARSALIB):
	$(MAKE) --directory=$(FARSADIR)/FaRSA/src

# Dependencies for libraries
$(FaRSALIB): $(wildcard $(FARSADIR)/FaRSA/src/*.hpp) $(wildcard $(FARSADIR)/FaRSA/src/*.cpp)

# Rule for objects
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Objects depend on headers and sources
$(objects): $(headers)
$(objects): $(sources)

# Clean
clean:
	rm -f $(objects) $(depends)

# Very clean
veryclean: clean
	rm -f $(EXES)
