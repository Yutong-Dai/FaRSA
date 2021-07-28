// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __FARSAITERATIONQUANTITIES_HPP__
#define __FARSAITERATIONQUANTITIES_HPP__

#include <ctime>
#include <memory>
#include <string>
#include <vector>
#include <iostream>

#include "FaRSAOptions.hpp"
#include "FaRSAPoint.hpp"
#include "FaRSAProblem.hpp"
#include "FaRSAReporter.hpp"
#include "FaRSAVector.hpp"

namespace FaRSA
{

/**
 * Forward declarations
 */
class Options;
class Point;
class Problem;
class Reporter;
class Vector;

/**
 * Quantities class
 */
class Quantities
{

public:
  /** @name Constructors */
  //@{
  /**
   * Declare Quantities
   */
  Quantities();
  //@}

  /** @name Destructor */
  //@{
  /**
   * Delete data
   */
  ~Quantities();
  //@}

  /** @name Options handling methods */
  //@{
  /**
   * Add options
   * \param[in,out] options is pointer to Options object from FaRSA
   * \param[in] reporter is pointer to Reporter object from FaRSA
   */
  void addOptions(Options* options,
                  const Reporter* reporter);
  /**
   * Set options
   * \param[in] options is pointer to Options object from FaRSA
   * \param[in] reporter is pointer to Reporter object from FaRSA
   */
  void getOptions(const Options* options,
                  const Reporter* reporter);
  //@}

  /** @name Initialization methods */
  //@{
  /**
   * Initialize quantities
   * \param[in] problem is pointer to Problem object
   */
  bool initialize(const std::shared_ptr<Problem> problem);
  //@}

  /** @name Get methods */
  //@{
  /**
   * Start time
   * \return start time that was set
   */
  inline clock_t const startTime() const { return start_time_; };
  /**
   * \return end time that was set
   */
  inline clock_t const endTime() const { return end_time_; };
  /**
   * Evaluation time
   * \return problem function evaluation time that was set
   */
  inline clock_t const evaluationTime() const { return evaluation_time_; };
  /**
   * Get CPU time limit
   * \return CPU time limit
   */
  inline double const cpuTimeLimit() const { return cpu_time_limit_; };
  /**
   * Get scaling threshold
   * \return scaling threshold
   */
  inline double const scalingThreshold() const { return scaling_threshold_; };
  /**
   * Get stepsize
   * \return current stepsize
   */
  inline double const stepsize() const { return stepsize_; };
  /**
   * Function evaluation counter
   * \return function evaluations performed so far
   */
  inline int const functionCounter() const { return function_counter_; };
  /**
   * Function evaluation limit
   * \return function evaluation limit
   */
  inline int const functionEvaluationLimit() const { return function_evaluation_limit_; };
  /**
   * Gradient evaluation counter
   * \return gradient evaluations performed so far
   */
  inline int const gradientCounter() const { return gradient_counter_; };
  /**
   * Gradient evaluation limit
   * \return gradient evaluation limit
   */
  inline int const gradientEvaluationLimit() const { return gradient_evaluation_limit_; };
  /**
   * Iteration counter
   * \return iterations performed so far
   */
  inline int const iterationCounter() const { return iteration_counter_; };
  /**
   * Get problem size
   * \return number of variables
   */
  inline int const numberOfVariables() const { return number_of_variables_; };
  /**
   * Get pointer to current iterate
   * \return pointer to Point representing current iterate
   */
  inline std::shared_ptr<Point> currentIterate() { return current_iterate_; };
  /**
   * Get pointer to trial iterate
   * \return pointer to Point representing trial iterate
   */
  inline std::shared_ptr<Point> trialIterate() { return trial_iterate_; };
  /**
   * Get direction
   * \return pointer to Vector representing search direction
   */
  inline std::shared_ptr<Vector> direction() { return direction_; };

  /**
   * @brief Get stationarity tolerance
   * \return double const stationarity_tolerance_
   */
  inline double const stationarityTolerance() const { return stationarity_tolerance_; };
  /**
   * @brief Get stationarity iteration limit
   * 
   * \return double const iteration_limit_
   */
  inline double const iterationLimit() const{ return iteration_limit_;} ;
  /**
   * @brief Get iterate norm tolerance
   * 
   * \return double const iterate_norm_tolerance_
   */
  inline double const iterateNormTolerance() const{ return iterate_norm_tolerance_;} ;

  inline double const linesearchArmijoEta() const{ return linesearch_armijo_eta_;} ;

  inline double const linesearchStepsizeDecreaseFactor() const{ return linesearch_stepsize_decrease_factor_;} ;

  inline int const linesearchMaxBacktrack() const{ return linesearch_max_backtrack_;} ;

  inline double const kappa1Max() const{ return kappa1_max_;} ;
  inline double const kappa1min() const{ return kappa1_min_;} ;
  inline double const kappa2Max() const{ return kappa2_max_;} ;
  inline double const kappa2min() const{ return kappa2_min_;} ;
  inline double const kappaIncreaseFactor() const{ return kappa_increase_factor_;} ;
  inline double const kappaDecreaseFactor() const{ return kappa_decrease_factor_;} ;



  //@}

  /** @name Set methods */
  //@{
  /**
   * Set current iterate pointer
   * \param[in] iterate is pointer to Point to represent current iterate
   */
  inline void setCurrentIterate(const std::shared_ptr<Point> iterate) { current_iterate_ = iterate; };
  /**
   * Set trial iterate pointer
   * \param[in] trial_iterate is pointer to Point to represent trial iterate
   */
  inline void setTrialIterate(const std::shared_ptr<Point> trial_iterate) { trial_iterate_ = trial_iterate; };
  /**
   * Set trial iterate pointer to current iterate pointer
   */
  inline void setTrialIterateToCurrentIterate() { trial_iterate_ = current_iterate_; };
  /**
   * Set stepsize
   * \param[in] stepsize is new value to represent stepsize
   */
  inline void setStepsize(double stepsize) { stepsize_ = stepsize; };
  //@}

  /** @name Increment methods */
  //@{
  /**
   * Increment evaluation time
   * \param[in] evaluation_time is amount to add to total problem function evaluation time
   */
  inline void incrementEvaluationTime(clock_t evaluation_time) { evaluation_time_ += evaluation_time; };
  /**
   * Increment function evaluation counter
   */
  inline void incrementFunctionCounter() { function_counter_++; };
  /**
   * Increment gradient evaluation counter
   */
  inline void incrementGradientCounter() { gradient_counter_++; };
  /**
   * Increment iteration counter
   */
  inline void incrementIterationCounter() { iteration_counter_++; };
  //@}

  /** @name Print methods */
  //@{
  /**
   * @brief print values of private members that are loaded from the Option instance.
   * 
   * \param[in] reporter is pointer to Reporter object from FaRSA
   */
  void print(const Reporter* reporter);
  /**
   * Get iteration header string
   * \return string of header values
   */
  std::string iterationHeader();
  /**
   * Get iteration null values string
   * \return string of null values
   */
  std::string iterationNullValues();
  //@}

  /** @name Print methods */
  //@{
  /**
   * Print header
   * \param[in] reporter is pointer to Reporter object from FaRSA
   */
  void printHeader(const Reporter* reporter);
  /**
   * Print iteration values
   * \param[in] reporter is pointer to Reporter object from FaRSA
   */
  void printIterationValues(const Reporter* reporter);
  /**
   * Print footer
   * \param[in] reporter is pointer to Reporter object from FaRSA
   */
  void printFooter(const Reporter* reporter);
  //@}

  /** @name Finalization method */
  //@{
  /**
   * Finalize
   */
  void finalize();
  //@}

private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Copy constructor
   */
  Quantities(const Quantities&);
  /**
   * Overloaded equals operator
   */
  void operator=(const Quantities&);
  //@}

  /** @name Private members */
  //@{
  clock_t end_time_;
  clock_t evaluation_time_;
  clock_t start_time_;
  double stepsize_;
  int function_counter_;
  int gradient_counter_;
  int iteration_counter_;
  int number_of_variables_;
  std::shared_ptr<Point> current_iterate_;
  std::shared_ptr<Point> trial_iterate_;
  std::shared_ptr<Vector> direction_;
  std::vector<int> groups_free_;
  std::vector<int> groups_zero_;
  //@}

  /** @name Private members (set by options) */
  //@{
  double scaling_threshold_;
  // for termination
  double cpu_time_limit_;
  double iterate_norm_tolerance_;
  double stationarity_tolerance_;
  int function_evaluation_limit_;
  int gradient_evaluation_limit_;
  int iteration_limit_;
  // for backtrack-linesearch
  double linesearch_armijo_eta_;
  double linesearch_stepsize_decrease_factor_;
  int linesearch_max_backtrack_;
  // for space partition using (FaRSAGroup)
  double kappa1_max_;
  double kappa1_min_;
  double kappa2_max_;
  double kappa2_min_;
  double kappa_increase_factor_;
  double kappa_decrease_factor_;
  //@}

}; // end Quantities

} // namespace FaRSA

#endif /* __FARSAITERATIONQUANTITIES_HPP__ */
