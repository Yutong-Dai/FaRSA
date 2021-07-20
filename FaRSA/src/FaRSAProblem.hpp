// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __FARSAPROBLEM_HPP__
#define __FARSAPROBLEM_HPP__

#include <vector>

namespace FaRSA
{

/**
 * Problem class
 */
class Problem
{

public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   */
  Problem(){};
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  virtual ~Problem(){};

  /** @name Get methods */
  //@{
  /**
   * Number of groups
   * \param[out] n_g is the number of groups, an integer (return value)
   * \return indicator of success (true) or failure (false)
   */
  inline bool numberOfGroups(int& n_g) { n_g = groups_.size(); return true; };
  /**
   * Number of variables
   * \param[out] n is the number of variables, an integer (return value)
   * \return indicator of success (true) or failure (false)
   */
  inline bool numberOfVariables(int& n) { n = number_of_variables_; return true; };
  /**
   * Returns initial point
   * \param[out] x is the initial point/iterate, a double array (return value)
   */
  virtual bool initialPoint(double* x) = 0;
  //@}

  /** @name Evaluate methods */
  //@{
  /**
   * Evaluates objective
   * \param[in] x is a given point/iterate, a constant double array
   * \param[out] f is the objective value at "x", a double (return value)
   */
  virtual bool evaluateObjective(const double* x,
                                 double& f) = 0;
  /**
   * Evaluates gradient
   * \param[in] x is a given point/iterate, a constant double array
   * \param[out] g is the gradient value at "x", a double array (return value)
   */
  virtual bool evaluateGradient(const double* x,
                                double* g) = 0;
  /**
   * Evaluates Hessian-vector product
   * \param[in] x is a given point/iterate, a constant double array
   * \param[in] groups is a vector of group indices
   * \param[in] v is a given vector, a constant double array
   * \param[out] Hv is the Hessian value at "x" times "v", a double array (return value)
   */
  virtual bool evaluateHessianVectorProduct(const double* x,
                                            const std::vector<int> groups,
                                            const double* v,
                                            double* Hv) = 0;
  //@}

  /** @name Finalize methods */
  //@{
  /**
   * Finalizes solution
   * \param[in] x is the final point/iterate, a constant double array
   * \param[in] f is the objective value at "x", a constant double
   * \param[in] g is the gradient value at "x", a constant double array
   */
  virtual bool finalizeSolution(const double* x,
                                double f,
                                const double* g) = 0;
  //@}

protected:
  /** @name Protected members */
  //@{
  int number_of_variables_;                /**< Number of variables */
  std::vector< std::vector<int> > groups_; /**< Group data          */
  //@}

private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Copy constructor
   */
  Problem(const Problem&);
  /**
   * Overloaded equals operator
   */
  void operator=(const Problem&);
  //@}

}; // end Problem

} // namespace FaRSA

#endif /* __FARSAPROBLEM_HPP__ */
