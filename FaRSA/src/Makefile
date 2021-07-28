# Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
#
# This code is published under the ??? License.
#
# Author(s) : Frank E. Curtis, Daniel P. Robinson

# C++ compiler
CXX = g++

# C++ compiler flags
CXXFLAGS = -g -Wall -std=c++11

# Library utility command
AR = ar rv

# Set sources, etc.
headers = $(wildcard *.hpp)
sources = $(wildcard *.cpp)
objects = $(sources:.cpp=.o)
depends = $(sources:.cpp=.d)
ofiles = $(wildcard *.o)
# Library
LIB = libFaRSA.a

# Rule for all
all: $(LIB)

# Create library
$(LIB): $(objects)
	$(AR) $@ $(objects)

# Rule for objects
.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Objects depend on headers and sources
$(objects): $(headers)
$(objects): $(sources)

# Clean
clean:
	rm -f $(objects) $(depends)

# Very clean
veryclean: clean
	rm -f $(LIB)

# for build efficiency
# single:
# 	@read -p "Filename : " Filename \
# 	&& $(CXX) $(CXXFLAGS) -c $${Filename}.cpp -o $${Filename}.o
single:
	$(CXX) $(CXXFLAGS) -c $(ARGS).cpp -o $(ARGS).o
rebuild:
	ar rv $(LIB) $(ofiles)