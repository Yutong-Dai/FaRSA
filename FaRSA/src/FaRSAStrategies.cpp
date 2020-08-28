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
                           "direction_computation",
                           "CuttingPlane",
                           "Direction computation strategy to use.\n"
                           "Default     : CuttingPlane.");
  options->addStringOption(reporter,
                           "line_search",
                           "WeakWolfe",
                           "Line search strategy to use.\n"
                           "Default     : WeakWolfe.");

  // Add options for direction computation strategies
  std::shared_ptr<DirectionComputation> direction_computation;
  direction_computation = std::make_shared<DirectionComputationProximalGradient>();
  direction_computation->addOptions(options, reporter);
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
  std::string direction_computation_name;
  std::string line_search_name;

  // Read integer options
  options->valueAsString(reporter, "direction_computation", direction_computation_name);
  options->valueAsString(reporter, "line_search", line_search_name);

  // Set direction computation strategy
  if (direction_computation_name.compare("ProximalGradient") == 0) {
    direction_computation_ = std::make_shared<DirectionComputationProximalGradient>();
  }
  else {
    direction_computation_ = std::make_shared<DirectionComputationProximalGradient>();
  }

  // Set line search strategy
  if (line_search_name.compare("Backtracking") == 0) {
    line_search_ = std::make_shared<LineSearchBacktracking>();
  }
  else {
    line_search_ = std::make_shared<LineSearchBacktracking>();
  }

  // Set direction computation options
  direction_computation_->getOptions(options, reporter);

  // Set line search options
  line_search_->getOptions(options, reporter);

} // end getOptions

// Initialize
void Strategies::initialize(const Options* options,
                            Quantities* quantities,
                            const Reporter* reporter)
{

  // Initialize direction computation
  direction_computation_->initialize(options, quantities, reporter);

  // Initialize line search
  line_search_->initialize(options, quantities, reporter);

} // end initialize

// Set iteration header
void Strategies::setIterationHeader()
{

  // Set iteration header string based on strategy objects
  iteration_header_ = "";
  if (direction_computation_->iterationHeader().length() > 0) {
    iteration_header_ += " ";
    iteration_header_ += direction_computation_->iterationHeader();
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
  reporter->printf(R_SOLVER, R_BASIC, "Direction computation strategy..... : %s\n"
                                      "Line search strategy............... : %s\n",
                   direction_computation_->name().c_str(),
                   line_search_->name().c_str());

} // end printHeader

// Print footer
void Strategies::printFooter(const Reporter* reporter) {}

} // namespace FaRSA
