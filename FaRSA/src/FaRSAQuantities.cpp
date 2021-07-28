// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#include <cmath>

#include "FaRSADefinitions.hpp"
#include "FaRSAQuantities.hpp"

namespace FaRSA
{

// Constructor
Quantities::Quantities()
  : evaluation_time_(0),
    stepsize_(0.0),
    function_counter_(0),
    gradient_counter_(0),
    iteration_counter_(0),
    number_of_variables_(0),
    scaling_threshold_(1.0),
    cpu_time_limit_(FARSA_DOUBLE_INFINITY),
    function_evaluation_limit_(1),
    gradient_evaluation_limit_(1)
{
  start_time_ = clock();
  end_time_ = start_time_;
  current_iterate_.reset();
  trial_iterate_.reset();
  direction_.reset();
}

// Destructor
Quantities::~Quantities() {}

// Add options
void Quantities::addOptions(Options* options,
                            const Reporter* reporter)
{

  // Add bool options

  // Add double options
  options->addDoubleOption(reporter,
                           "cpu_time_limit",
                           1e+04,
                           0.0,
                           FARSA_DOUBLE_INFINITY,
                           "Limit on the number of CPU seconds.  This limit is only checked\n"
                           "              at the beginning of an iteration, so the true CPU time limit\n"
                           "              also depends on the time required to a complete an iteration.\n"
                           "Default     : 1e+04.");
  options->addDoubleOption(reporter,
                           "scaling_threshold",
                           1e+02,
                           0.0,
                           FARSA_DOUBLE_INFINITY,
                           "Threshold for determining objective scaling.  If norm of gradient\n"
                           "              at the initial point is greater than this value, then the objective\n"
                           "              is scaled so that the initial gradient norm is at this value.\n"
                           "Default     : 1e+02.");
  options->addDoubleOption(reporter,
                           "iterate_norm_tolerance",
                           1e+20,
                           0.0,
                           FARSA_DOUBLE_INFINITY,
                           "Tolerance for determining divergence of the algorithm iterates.\n"
                           "              If the norm of an iterate is larger than this tolerance times\n"
                           "              the maximum of 1.0 and the norm of the initial iterate, then\n"
                           "              the algorithm terminates with a message of divergence.\n"
                           "Default     : 1e+20.");
  options->addDoubleOption(reporter,
                           "stationarity_tolerance",
                           1e-04,
                           0.0,
                           FARSA_DOUBLE_INFINITY,
                           "Tolerance for determining stationarity.  If the stationarity\n"
                           "              measure falls below this tolerance, then the algorithm\n"
                           "              terminates with a message of stationarity.\n"
                           "Default     : 1e-04.");

  // Add integer options
  options->addIntegerOption(reporter,
                            "function_evaluation_limit",
                            1e+05,
                            0,
                            FARSA_INT_INFINITY,
                            "Limit on the number of function evaluations performed.\n"
                            "Default     : 1e+05.");
  options->addIntegerOption(reporter,
                            "gradient_evaluation_limit",
                            1e+05,
                            0,
                            FARSA_INT_INFINITY,
                            "Limit on the number of gradient evaluations performed.\n"
                            "Default     : 1e+05.");
  
  options->addIntegerOption(reporter,
                            "iteration_limit",
                            1e+04,
                            0,
                            FARSA_INT_INFINITY,
                            "Limit on the number of iterations that will be performed.\n"
                            "              Note that each iteration might involve inner iterations.\n"
                            "Default     : 1e+04.");

} // end addOptions

// Set options
void Quantities::getOptions(const Options* options,
                            const Reporter* reporter)
{
  // set bool options
  
  // Set double options
  options->valueAsDouble(reporter, "iterate_norm_tolerance", iterate_norm_tolerance_);
  options->valueAsDouble(reporter, "stationarity_tolerance", stationarity_tolerance_);
  options->valueAsDouble(reporter, "cpu_time_limit", cpu_time_limit_);
  options->valueAsDouble(reporter, "scaling_threshold", scaling_threshold_);

  // Set integer options
  options->valueAsInteger(reporter, "function_evaluation_limit", function_evaluation_limit_);
  options->valueAsInteger(reporter, "gradient_evaluation_limit", gradient_evaluation_limit_);
  options->valueAsInteger(reporter, "iteration_limit", iteration_limit_);

} // end getOptions

// Initialization
bool Quantities::initialize(const std::shared_ptr<Problem> problem)
{

  // Start clock
  start_time_ = clock();
  end_time_ = start_time_;

  // Initialize counters
  evaluation_time_ = 0;
  function_counter_ = 0;
  gradient_counter_ = 0;
  iteration_counter_ = 0;

  // Declare success boolean
  bool success = true;

  // Declare integer
  int n;

  // Get number of variables
  success = problem->numberOfVariables(n);

  // Check for success
  if (!success) {
    return false;
  }

  // Set number of variables
  number_of_variables_ = n;

  // Declare vector
  std::shared_ptr<Vector> v(new Vector(number_of_variables_));

  // Get initial point
  success = problem->initialPoint(v->valuesModifiable());

  // Check for success
  if (!success) {
    return false;
  }

  // Declare iterate
  std::shared_ptr<Point> initial_iterate(new Point(problem, v, 1.0));

  // Set initial point
  current_iterate_ = initial_iterate;

  // Initialize direction
  direction_ = std::make_shared<Vector>(number_of_variables_);

  // Initialize stepsize
  stepsize_ = 0.0;

  // Return
  return success;

} // end initialize

// Iteration header string
std::string Quantities::iterationHeader()
{
  return "  Iter.  Objective ";
}

// Iteration null values
std::string Quantities::iterationNullValues()
{
  return " ------ -----------";
}

// Print header
void Quantities::printHeader(const Reporter* reporter)
{

  // Print header
  reporter->printf(R_SOLVER, R_BASIC, "Number of variables................ : %d\n", number_of_variables_);

} // end printHeader

// Print iteration values
void Quantities::printIterationValues(const Reporter* reporter)
{

  // Print iteration values
  reporter->printf(R_SOLVER, R_PER_ITERATION, " %6d %+.4e", iteration_counter_, current_iterate_->objective());

} // end printIterationValues

// Print footer
void Quantities::printFooter(const Reporter* reporter)
{

  // Print quantities footer
  reporter->printf(R_SOLVER, R_BASIC, "\n\n"
                                  "Objective.......................... : %e\n"
                                  "Objective (unscaled)............... : %e\n"
                                  "\n"
                                  "Number of iterations............... : %d\n"
                                  "Number of function evaluations..... : %d\n"
                                  "Number of gradient evaluations..... : %d\n"
                                  "\n"
                                  "CPU seconds........................ : %f\n"
                                  "CPU seconds in FaRSA............... : %f\n"
                                  "CPU seconds in evaluations......... : %f\n",
                   current_iterate_->objective(),
                   current_iterate_->objectiveUnscaled(),
                   iteration_counter_,
                   function_counter_,
                   gradient_counter_,
                   (end_time_ - start_time_) / (double)CLOCKS_PER_SEC,
                   (end_time_ - start_time_ - evaluation_time_) / (double)CLOCKS_PER_SEC,
                   evaluation_time_ / (double)CLOCKS_PER_SEC);

} // end printFooter

// Finalization
void Quantities::finalize()
{

  // Set end time
  end_time_ = clock();

} // end finalize

} // namespace FaRSA
