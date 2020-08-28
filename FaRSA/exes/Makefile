# Copyright (C) 2019 Frank E. Curtis
#
# This code is published under the MIT License.
#
# Author(s) : Frank E. Curtis

# C++ compiler
CXX = g++

# C++ compiler flags
CXXFLAGS = -g -Wall -std=c++11

# Set sources, etc.
sources = $(wildcard *.cpp)
objects = $(sources:.cpp=.o)
depends = $(sources:.cpp=.d)

# Executable(s)
EXES = $(sources:.cpp=)

# Libraries
FaRSALIB         = $(FARSADIR)/FaRSA/src/libFaRSA.a
FaRSAProblemsLIB = $(FARSADIR)/FaRSA/problems/libFaRSAProblems.a

# Includes
INCLUDES = -I $(FARSADIR)/FaRSA/src -I $(FARSADIR)/FaRSA/problems

# Rule for all
all: $(EXES)

# Dependence for executable
$(EXES): % : %.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(FaRSALIB) $(FaRSAProblemsLIB) -L $(LAPACKDIR) -ldl -lblas -llapack

# Dependencies for executable
$(EXES): $(FaRSALIB) $(FaRSAProblemsLIB)

# Rules for libraries
$(FaRSALIB):
	$(MAKE) --directory=$(FARSADIR)/FaRSA/src
$(FaRSAProblemsLIB):
	$(MAKE) --directory=$(FARSADIR)/FaRSA/problems

# Dependencies for libraries
$(FaRSALIB): $(wildcard $(FARSADIR)/FaRSA/src/*.hpp) $(wildcard $(FARSADIR)/FaRSA/src/*.cpp)
$(FaRSAProblemsLIB): $(wildcard $(FARSADIR)/FaRSA/problems/*.hpp) $(wildcard $(FARSADIR)/FaRSA/problems/*.cpp)

# Rule for .cpp.o
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Objects depend on sources
$(objects): $(sources)

# Clean
clean:
	rm -f $(objects) $(depends)

# Very clean
veryclean: clean
	rm -f $(EXES)
