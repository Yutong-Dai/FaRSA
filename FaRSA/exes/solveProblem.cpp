// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#include <cstdio>
#include <cstring>
#include <string>

#include "FaRSAProblem.hpp"
#include "FaRSASolver.hpp"
#include "SimpleQuadratic.hpp"

using namespace FaRSA;

// Main function
int main(int argc, char* argv[])
{

  // Set usage string
  std::string usage("Usage: ./solveProblem ProblemName\n"
                    "       where ProblemName is name of problem in problems subdirectory.\n");

  // Check number of input arguments
  if (argc < 2) {
    printf("Too few arguments. Quitting.\n");
    printf("%s", usage.c_str());
    return 1;
  }

  // Declare problem dimension
  int const dimension = 10;

  // Declare problem
  std::shared_ptr<Problem> problem;
  if (strcmp(argv[1], "SimpleQuadratic") == 0) {
    problem = std::make_shared<SimpleQuadratic>(dimension);
  }
  else {
    printf("Invalid problem name. Quitting.\n");
    printf("%s", usage.c_str());
    return 1;
  }

  // Declare solver object
  FaRSASolver farsa;

  // Modify options from file
  farsa.options()->modifyOptionsFromFile(farsa.reporter(), "farsa.opt");

  // Optimize
  farsa.optimize(problem);

  // Return
  return 0;

} // end main
