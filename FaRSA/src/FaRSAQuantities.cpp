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

  options->addDoubleOption(reporter,
                           "linesearch_armijo_eta",
                           1e-04,
                           0.0,
                           1.0,
                           "Constant used in the Armijo backtrack-linesearch.\n"
                           "              Specifically, it reduce the directional derivative\n"
                           "              by a factor of eta: \eta * <gradient, d>.\n"
                           "Default     : 1e-04.");      
  options->addDoubleOption(reporter,
                           "linesearch_stepsize_decrease_factor",
                           0.8,
                           0.0,
                           1.0,
                           "Constant used to decrease the stepsize during the linesearch.\n"
                           "Default     : 0.8.");                                                
  options->addDoubleOption(reporter,
                           "kappa1_max",
                           1e+06,
                           0.0,
                           FARSA_DOUBLE_INFINITY,
                           "kappa1_max.\n"
                           "Default     : 1e+06.");
  options->addDoubleOption(reporter,
                           "kappa1_min",
                           1e-05,
                           0.0,
                           FARSA_DOUBLE_INFINITY,
                           "kappa1_min.\n"
                           "Default     : 1e-05."); 
  options->addDoubleOption(reporter,
                           "kappa2_max",
                           1e+05,
                           0.0,
                           FARSA_DOUBLE_INFINITY,
                           "kappa2_max.\n"
                           "Default     : 1e+05.");
  options->addDoubleOption(reporter,
                           "kappa2_min",
                           1e-06,
                           0.0,
                           FARSA_DOUBLE_INFINITY,
                           "kappa2_min.\n"
                           "Default     : 1e-06.");       
  options->addDoubleOption(reporter,
                            "kappa_increase_factor",
                            10.0,
                            1.0,
                            FARSA_DOUBLE_INFINITY,
                            "Factor to increase the kappa1 and kappa2.\n"
                            "Default     : 10.");
  options->addDoubleOption(reporter,
                            "kappa_decrease_factor",
                            0.1,
                            0.0,
                            1.0,
                            "Factor to decrease the kappa1 and kappa2.\n"
                            "Default     : 0.1.");                               
                                                                    

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
  options->addIntegerOption(reporter,
                            "linesearch_max_backtrack",
                            100,
                            0,
                            FARSA_INT_INFINITY,
                            "Limit on the number of backtrack-linesearch will be performed.\n"
                            "Default     : 100.");                                                                              

} // end addOptions

// Set options
void Quantities::getOptions(const Options* options,
                            const Reporter* reporter)
{
  // set bool options
  
  // Set double options
  options->valueAsDouble(reporter, "cpu_time_limit", cpu_time_limit_);
  options->valueAsDouble(reporter, "iterate_norm_tolerance", iterate_norm_tolerance_);
  options->valueAsDouble(reporter, "stationarity_tolerance", stationarity_tolerance_);
  options->valueAsDouble(reporter, "scaling_threshold", scaling_threshold_);
  options->valueAsDouble(reporter, "linesearch_armijo_eta", linesearch_armijo_eta_);
  options->valueAsDouble(reporter, "linesearch_stepsize_decrease_factor", linesearch_stepsize_decrease_factor_);
  options->valueAsDouble(reporter, "kappa1_max", kappa1_max_);
  options->valueAsDouble(reporter, "kappa1_min", kappa1_min_);
  options->valueAsDouble(reporter, "kappa2_max", kappa2_max_);
  options->valueAsDouble(reporter, "kappa2_min", kappa2_min_);
  options->valueAsDouble(reporter, "kappa_increase_factor", kappa_increase_factor_);
  options->valueAsDouble(reporter, "kappa_decrease_factor", kappa_decrease_factor_);          

  // Set integer options
  options->valueAsInteger(reporter, "function_evaluation_limit", function_evaluation_limit_);
  options->valueAsInteger(reporter, "gradient_evaluation_limit", gradient_evaluation_limit_);
  options->valueAsInteger(reporter, "iteration_limit", iteration_limit_);
  options->valueAsInteger(reporter, "linesearch_max_backtrack", linesearch_max_backtrack_);
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

// Print all values of private memebers set by the Option class
void Quantities::print(const Reporter* reporter){
  reporter->printf(R_SOLVER, R_BASIC, "---------------------- Termination Conditions -----------------\n");
  reporter->printf(R_SOLVER, R_BASIC, "cpu time limit (seconds)...................... : %+.4e\n", cpu_time_limit_);
  reporter->printf(R_SOLVER, R_BASIC, "iterate norm tolerance........................ : %+.4e\n", iterate_norm_tolerance_);
  reporter->printf(R_SOLVER, R_BASIC, "stationarity tolerance........................ : %+.4e\n", stationarity_tolerance_);
  reporter->printf(R_SOLVER, R_BASIC, "function evaluation limit_.................... : %d\n", function_evaluation_limit_);
  reporter->printf(R_SOLVER, R_BASIC, "gradient evaluation limit..................... : %d\n", gradient_evaluation_limit_);
  reporter->printf(R_SOLVER, R_BASIC, "iteration limit............................... : %d\n", iteration_limit_);
  reporter->printf(R_SOLVER, R_BASIC, "---------------------- Algorithmic Choices ---------------------\n");
  reporter->printf(R_SOLVER, R_BASIC, "scaling threshold............................. : %+.4e\n", scaling_threshold_);
  reporter->printf(R_SOLVER, R_BASIC, "linesearch armijo eta......................... : %+.4e\n", linesearch_armijo_eta_);
  reporter->printf(R_SOLVER, R_BASIC, "linesearch stepsize decrease factor........... : %+.4e\n", linesearch_stepsize_decrease_factor_);
  reporter->printf(R_SOLVER, R_BASIC, "linesearch max backtrack...................... : %d\n", linesearch_max_backtrack_);
  reporter->printf(R_SOLVER, R_BASIC, "kappa1 max.................................... : %+.4e\n", kappa1_max_);
  reporter->printf(R_SOLVER, R_BASIC, "kappa1 min.................................... : %+.4e\n", kappa1_min_);
  reporter->printf(R_SOLVER, R_BASIC, "kappa2 max.................................... : %+.4e\n", kappa2_max_);
  reporter->printf(R_SOLVER, R_BASIC, "kappa2 min.................................... : %+.4e\n", kappa2_min_);
  reporter->printf(R_SOLVER, R_BASIC, "kappa increase factor......................... : %+.4e\n", kappa_increase_factor_);
  reporter->printf(R_SOLVER, R_BASIC, "kappa decrease factor......................... : %+.4e\n", kappa_decrease_factor_);

}// end print

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
