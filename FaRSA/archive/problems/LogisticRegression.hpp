// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

// Description : Implementation for FaRSA of the objective
//                 f(x) = logistic regression + group regularizer

#ifndef __LOGISTICREGRESSION_HPP__
#define __LOGISTICREGRESSION_HPP__

#include <vector>

#include "FaRSAMatrix.hpp"
#include "FaRSAProblem.hpp"
#include "FaRSAVector.hpp"

using namespace FaRSA;

/**
 * LogisticRegression class
 */
class LogisticRegression : public Problem
{

public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   */
  LogisticRegression(char* features_file,
                     char* labels_file,
                     char* groups_file,
                     char* initial_point_file);
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  ~LogisticRegression();
  //@}

  /** @name Get methods */
  //@{
  /**
   * Number of data points
   * \param[out] N is the number of data points, an integer (return value)
   * \return indicator of success (true) or failure (false)
   */
  inline bool numberOfDataPoints(int& N) { N = number_of_data_points_; return true; };
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
  LogisticRegression(const LogisticRegression&);
  /**
   * Overloaded equals operator
   */
  void operator=(const LogisticRegression&);
  //@}

  /** @name Private members */
  //@{
  int number_of_data_points_; /**< Number of data points */
  Vector initial_point_;      /**< Initial point         */
  Matrix features_;           /**< Feature data          */
  Vector labels_;             /**< Label data            */
  //@}

  /** @name Private members */
  //@{
  void setGroupsFromFile(char* groups_file);
  //@}

}; // end LogisticRegression

#endif /* __LOGISTICREGRESSION_HPP__ */
