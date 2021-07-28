// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __FARSASOLVER_HPP__
#define __FARSASOLVER_HPP__

#include <ctime>
#include <memory>

#include "FaRSAEnumerations.hpp"
#include "FaRSAOptions.hpp"
#include "FaRSAProblem.hpp"
#include "FaRSAQuantities.hpp"
#include "FaRSAReporter.hpp"
#include "FaRSAStrategies.hpp"

namespace FaRSA
{

/**
 * Forward declarations
 */
class Options;
class Problem;
class Quantities;
class Reporter;
class Strategies;

/**
 * FaRSASolver class
 */
class FaRSASolver
{

public:
  /** @name Constructors */
  //@{
  /**
   * Construct FaRSASolver
   */
  FaRSASolver();
  //@}

  /** @name Destructor */
  //@{
  /**
   * Delete FaRSASolver
   */
  ~FaRSASolver();

  /** @name Get methods */
  //@{
  /**
   * Get function evaluation counter
   * \return function evaluations so far
   */
  inline int const functionEvaluations() const { return quantities_.functionCounter(); };
  /**
   * Get gradient evaluation counter
   * \return gradient evaluations so far
   */
  inline int const gradientEvaluations() const { return quantities_.gradientCounter(); };
  /**
   * Get iteration counter
   * \return iterations performed so far
   */
  inline int const iterations() const { return quantities_.iterationCounter(); };
  /**
   * Get number of variables
   * \return number of variables
   */
  inline int const numberOfVariables() const { return quantities_.numberOfVariables(); };
  /**
   * Get objective value
   * \return objective value of current iterate
   */
  inline double const objective() { return quantities_.currentIterate()->objectiveUnscaled(); };
  /**
   * Get time in evaluations
   * \return seconds between start and end time
   */
  inline double const time() const { return (quantities_.endTime() - quantities_.startTime()) / (double)CLOCKS_PER_SEC; };
  /**
   * Get time in evaluations
   * \return seconds to perform problem function evaluations
   */
  inline double const timeEvaluations() const { return quantities_.evaluationTime() / (double)CLOCKS_PER_SEC; };
  /**
   * Get time in FaRSA
   * \return seconds between start and end time not including problem function evaluation time
   */
  inline double const timeFaRSA() const { return (quantities_.endTime() - quantities_.startTime() - quantities_.evaluationTime()) / (double)CLOCKS_PER_SEC; };
  /**
   * Get status
   * \return current status of algorithm
   */
  inline FaRSA_Status const status() const { return status_; };
  /**
   * Get options
   * \return pointer to Options object
   */
  inline Options* options() { return &options_; };
  /**
   * Get reporter
   * \return pointer to Reporter object
   */
  inline Reporter* reporter() { return &reporter_; };
  /**
   * Get solution
   * \param[out] vector is the current iterate
   */
  void solution(double vector[]);
  //@}

  /** @name Set method */
  //@{
  /**
   * Set status
   * \param[in] status is new status to be set
   */
  inline void setStatus(FaRSA_Status status) { status_ = status; };
  //@}

  /** @name Optimize method */
  //@{
  /**
   * Optimize
   * \param[in] problem is a pointer to a Problem object
   */
  void optimize(const std::shared_ptr<Problem> problem);
  //@}

private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Copy constructor
   */
  FaRSASolver(const FaRSASolver&);
  /**
   * Overloaded equals operator
   */
  void operator=(const FaRSASolver&);
  //@}

  /** @name Private members */
  //@{
  FaRSA_Status status_;
  //@}

  /** @name Private members, objects */
  //@{
  Options options_;
  Quantities quantities_;
  Reporter reporter_;
  Strategies strategies_;
  //@}

  /** @name Private methods */
  //@{
  void addOptions();
  void evaluateFunctionsAtCurrentIterate();
  void printFooter();
  void printHeader();
  void printIterationHeader();
  void getOptions();
  //@}

}; // end FaRSASolver

} // namespace FaRSA

#endif /* __FARSASOLVER_HPP__ */
