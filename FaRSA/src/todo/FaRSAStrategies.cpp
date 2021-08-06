// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#include "FaRSAStrategies.hpp"
#include "FaRSADirectionComputationProximalGradient.hpp"
#include "FaRSALineSearchBacktracking.hpp"

namespace FaRSA
{

// Add options
void Strategies::addOptions(Options* options,
                            const Reporter* reporter)
{

  // Add string options
  options->addStringOption(reporter,
                           "direction_computation_first_order",
                           "ProximalGradient",
                           "Direction strategy used to optimize the first order variables.\n"
                           "Default     : ProximalGradient.");
  options->addStringOption(reporter,
                           "direction_computation_second_order",
                           "TruncatedNewton",
                           "Direction strategy used to optimize the second order variables.\n"
                           "Default     : TruncatedNewton.");  
  options->addStringOption(reporter,
                           "line_search",
                           "Backtracking",
                           "Line search strategy to use.\n"
                           "Default     : Backtracking.");

  // Add options for space partition strategies

  // ADD NEW DIRECTION COMPUTATION STRATEGIES HERE AND IN SWITCH BELOW //


  // Add options for direction computation strategies
  std::shared_ptr<DirectionComputation> direction_computation_first_order;
  direction_computation_first_order = std::make_shared<DirectionComputationProximalGradient>();
  direction_computation_first_order->addOptions(options, reporter);
  // ADD NEW DIRECTION COMPUTATION STRATEGIES HERE AND IN SWITCH BELOW //

  // Add options for line search strategies
  std::shared_ptr<LineSearch> line_search;
  line_search = std::make_shared<LineSearchBacktracking>();
  line_search->addOptions(options, reporter);
  // ADD NEW LINE SEARCH STRATEGIES HERE AND IN SWITCH BELOW //

} // end addOptions

// Set options
void Strategies::getOptions(const Options* options,
                            const Reporter* reporter)
{

  // Declare strategy names
  std::string direction_computation_first_order_name;
  std::string line_search_name;

  // Read integer options
  options->valueAsString(reporter, "direction_computation_first_order", direction_computation_first_order_name);
  options->valueAsString(reporter, "line_search", line_search_name);

  // Set direction computation strategy
  if (direction_computation_first_order_name.compare("ProximalGradient") == 0) {
    direction_computation_first_order_ = std::make_shared<DirectionComputationProximalGradient>();
  }
  else {
    direction_computation_first_order_ = std::make_shared<DirectionComputationProximalGradient>();
  }

  // Set line search strategy
  if (line_search_name.compare("Backtracking") == 0) {
    line_search_ = std::make_shared<LineSearchBacktracking>();
  }
  else {
    line_search_ = std::make_shared<LineSearchBacktracking>();
  }

  // Set direction computation options
  direction_computation_first_order_->getOptions(options, reporter);

  // Set line search options
  line_search_->getOptions(options, reporter);

} // end getOptions

// Initialize
void Strategies::initialize(const Options* options,
                            Quantities* quantities,
                            const Reporter* reporter)
{

  // Initialize direction computation
  direction_computation_first_order_->initialize(options, quantities, reporter);

  // Initialize line search
  line_search_->initialize(options, quantities, reporter);

} // end initialize

// Set iteration header
void Strategies::setIterationHeader()
{

  // Set iteration header string based on strategy objects
  iteration_header_ = "";
  if (direction_computation_first_order_->iterationHeader().length() > 0) {
    iteration_header_ += " ";
    iteration_header_ += direction_computation_first_order_->iterationHeader();
  } // end if
  if (line_search_->iterationHeader().length() > 0) {
    iteration_header_ += " ";
    iteration_header_ += line_search_->iterationHeader();
  } // end if

} // end setIterationHeader

// Print header
void Strategies::printHeader(const Reporter* reporter)
{

  // Print header
  reporter->printf(R_SOLVER, R_BASIC, "Direction computation strategy (1st order)..... : %s\n"
                                      "Line search strategy........................... : %s\n",
                   direction_computation_first_order_->name().c_str(),
                   line_search_->name().c_str());

} // end printHeader

// Print footer
void Strategies::printFooter(const Reporter* reporter) {}

} // namespace FaRSA
