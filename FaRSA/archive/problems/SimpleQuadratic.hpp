// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

// Description : Implementation for FaRSA of the objective
//                 f(x) = sum_{i=1..n} i*x_i^2
//               with initial point
//                 x_i = 1 for all i = 1..n
// Notes       : Optimal value: 0.0

#ifndef __SIMPLEQUADRATIC_HPP__
#define __SIMPLEQUADRATIC_HPP__

#include <vector>

#include "FaRSAProblem.hpp"

using namespace FaRSA;

/**
 * SimpleQuadratic class
 */
class SimpleQuadratic : public Problem
{

public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   */
  SimpleQuadratic(int n);
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  ~SimpleQuadratic();
  //@}

  /** @name Get methods */
  //@{
  /**
   * Initial point
   * \param[out] x is the initial point/iterate, a double array (return value)
   * \return indicator of success (true) or failure (false)
   */
  bool initialPoint(double* x);
  //@}

  /** @name Evaluate methods */
  //@{
  /**
   * Evaluates objective
   * \param[in] x is a given point/iterate, a constant double array
   * \param[out] f is the objective value at "x", a double (return value)
   * \return indicator of success (true) or failure (false)
   */
  bool evaluateObjective(const double* x,
                         double& f);
  /**
   * Evaluates gradient
   * \param[in] x is a given point/iterate, a constant double array
   * \param[out] g is the gradient value at "x", a double array (return value)
   * \return indicator of success (true) or failure (false)
   */
  bool evaluateGradient(const double* x,
                        double* g);
  /**
   * Evaluates Hessian-vector product
   * \param[in] x is a given point/iterate, a constant double array
   * \param[in] groups is a vector of group indices
   * \param[in] v is a given vector, a constant double array
   * \param[out] Hv is the product of the Hessian and "v", a double array (return value)
   * \return indicator of success (true) or failure (false)
   */
  bool evaluateHessianVectorProduct(const double* x,
                                    const std::vector<int> groups,
                                    const double* v,
                                    double* Hv);
  //@}

  /** @name Finalize methods */
  //@{
  /**
   * Finalizes solution
   * \param[in] x is the final point/iterate, a constant double array
   * \param[in] f is the objective value at "x", a constant double
   * \param[in] g is the gradient value at "x", a constant double array
   * \return indicator of success (true) or failure (false)
   */
  bool finalizeSolution(const double* x,
                        double f,
                        const double* g);
  //@}

private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Copy constructor
   */
  SimpleQuadratic(const SimpleQuadratic&);
  /**
   * Overloaded equals operator
   */
  void operator=(const SimpleQuadratic&);
  //@}

}; // end SimpleQuadratic

#endif /* __SIMPLEQUADRATIC_HPP__ */
