// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#include <cmath>

#include "SimpleQuadratic.hpp"

// Constructor
SimpleQuadratic::SimpleQuadratic(int n)
{

  // Set number of variables
  number_of_variables_ = n;

  // Set groups
  for (int i = 0; i < number_of_variables_; i++) {
    groups_[i].push_back(i);
  }

} // end constructor

// Destructor
SimpleQuadratic::~SimpleQuadratic() {}

// Initial point
bool SimpleQuadratic::initialPoint(double* x)
{

  // Set initial point
  for (int i = 0; i < number_of_variables_; i++) {
    x[i] = 1.0;
  }

  // Return
  return true;

} // end initialPoint

// Objective value
bool SimpleQuadratic::evaluateObjective(const double* x,
                                        double& f)
{

  // Evaluate function
  f = 0.0;
  for (int i = 0; i < number_of_variables_; i++) {
    f += (double)(i+1) * pow(x[i],2.0);
  }

  // Return
  return true;

} // end evaluateObjective

// Gradient value
bool SimpleQuadratic::evaluateGradient(const double* x,
                                       double* g)
{

  // Evaluate gradient
  for (int i = 0; i < number_of_variables_; i++) {
    g[i] = (double)(i+1) * 2.0 * x[i];
  } // end for

  // Return
  return true;

} // end evaluateGradient

// Hessian-vector product
bool SimpleQuadratic::evaluateHessianVectorProduct(const double* x,
                                                   const std::vector<int> groups,
                                                   const double* v,
                                                   double* Hv)
{

  // Evaluate product
  for (int i = 0; i < number_of_variables_; i++) {
    Hv[i] = (double)(i+1) * 2.0;
  }

  // Return
  return true;

} // end evaluateHessianVectorProduct

// Finalize solution
bool SimpleQuadratic::finalizeSolution(const double* x,
                                       double f,
                                       const double* g)
{
  return true;
}
