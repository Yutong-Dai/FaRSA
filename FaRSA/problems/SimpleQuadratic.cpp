// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#include <cmath>

#include "SimpleQuadratic.hpp"

// Constructor
SimpleQuadratic::SimpleQuadratic(int n)
  : number_of_variables_(n) {}

// Destructor
SimpleQuadratic::~SimpleQuadratic() {}

// Number of variables
bool SimpleQuadratic::numberOfVariables(int& n)
{

  // Set number of variables
  n = number_of_variables_;

  // Return
  return true;

} // end numberOfVariables

// Initial point
bool SimpleQuadratic::initialPoint(int n,
                                   double* x)
{

  // Set initial point
  for (int i = 0; i < n; i++) {
    x[i] = 1.0;
  }

  // Return
  return true;

} // end initialPoint

// Objective value
bool SimpleQuadratic::evaluateObjective(int n,
                                        const double* x,
                                        double& f)
{

  // Evaluate function
  f = 0.0;
  for (int i = 0; i < n; i++) {
    f += (double)(i+1) * pow(x[i],2.0);
  }

  // Return
  return true;

} // end evaluateObjective

// Gradient value
bool SimpleQuadratic::evaluateGradient(int n,
                                       const double* x,
                                       double* g)
{

  // Evaluate gradient
  for (int i = 0; i < n; i++) {
    g[i] = (double)(i+1) * 2.0 * x[i];
  } // end for

  // Return
  return true;

} // end evaluateGradient

// Hessian-vector product
bool SimpleQuadratic::evaluateHessianVectorProduct(int n,
                                                   const double* x,
                                                   const double* v,
                                                   double* Hv)
{

  // Evaluate product
  for (int i = 0; i < n; i++) {
    Hv[i] = (double)(i+1) * 2.0;
  }

  // Return
  return true;

} // end evaluateHessianVectorProduct

// Finalize solution
bool SimpleQuadratic::finalizeSolution(int n,
                                       const double* x,
                                       double f,
                                       const double* g)
{
  return true;
}
