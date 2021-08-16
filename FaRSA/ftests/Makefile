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

FaRSAHeaders = $(wildcard ../src/*.hpp)
FaRSASource= $(wildcard ../src/*.cpp)

# Executable(s)
EXES = $(sources:.cpp=)

# Libraries
FaRSALIB = $(FARSADIR)/FaRSA/src/libFaRSA.a
FaRSAFUNCTIONLIB = $(FARSADIR)/FaRSA/functions/libFaRSAFunctions.a

# Includes
INCLUDES = -I $(FARSADIR)/FaRSA/src -I $(FARSADIR)/FaRSA/functions

# Rule for all
all: $(EXES)

# Rule for executable
$(EXES): % : %.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(FaRSALIB) $(FaRSAFUNCTIONLIB) -L $(LAPACKDIR) -ldl -lblas -llapack

# Dependencies for executable
$(EXES): $(FARSALIB) $(FaRSAFUNCTIONLIB)

# Rules for libraries
$(FARSALIB):
	$(MAKE) --directory=$(FARSADIR)/FaRSA/src
$(FaRSAFUNCTIONLIB):	
	$(MAKE) --directory=$(FARSADIR)/FaRSA/functions
# Dependencies for libraries
$(FaRSALIB): $(wildcard $(FARSADIR)/FaRSA/src/*.hpp) $(wildcard $(FARSADIR)/FaRSA/src/*.cpp)
$(FaRSAFUNCTIONLIB): $(wildcard $(FARSADIR)/FaRSA/functions/*.hpp) $(wildcard $(FARSADIR)/FaRSA/functions/*.cpp)

# Rule for objects
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Objects depend on headers and sources
$(objects): $(headers) $(FaRSAHeaders)
$(objects): $(sources) $(FaRSASource)

# Clean
clean:
	rm -f $(objects) $(depends)

# Very clean
veryclean: clean
	rm -f $(EXES)

# rebuild single test
single:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $(ARGS).cpp -o $(ARGS).o
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(ARGS) $(ARGS).o $(FaRSALIB) $(FaRSAFUNCTIONLIB) -L $(LAPACKDIR) -ldl -lblas -llapack