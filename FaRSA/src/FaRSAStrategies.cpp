// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#include "FaRSAStrategies.hpp"

#include "FaRSADirectionComputationProximalGradient.hpp"
#include "FaRSALineSearchBacktracking.hpp"
#include "FaRSAParameterUpdatePGStepsize.hpp"
#include "FaRSASpacePartitionFirstOrder.hpp"

namespace FaRSA
{
// Add options
void Strategies::addOptions(Options* options, const Reporter* reporter)
{
    // Add string options
    options->addStringOption(reporter, "space_partition", "FirstOrderPartition",
                             "Space partition strategy used to split groups into first order and "
                             "second order variable sets respectively.\n"
                             "Default     : FirstOrderPartition.");
    options->addStringOption(reporter, "direction_computation_first_order", "ProximalGradient",
                             "Direction strategy used to optimize the first order variables.\n"
                             "Default     : ProximalGradient.");
    options->addStringOption(reporter, "direction_computation_second_order", "TruncatedNewton",
                             "Direction strategy used to optimize the second order variables.\n"
                             "Default     : TruncatedNewton.");
    options->addStringOption(reporter, "line_search", "Backtracking",
                             "Line search strategy to use.\n"
                             "Default     : Backtracking.");

    // Add options for space partition strategies
    std::shared_ptr<SpacePartition> space_partition;
    space_partition = std::make_shared<SpacePartitionFirstOrder>();
    space_partition->addOptions(options, reporter);
    // ADD NEW SPACE PARTITION STRATEGIES HERE AND IN SWITCH BELOW //

    // Add options for direction computation strategies
    std::shared_ptr<DirectionComputation> direction_computation_first_order;
    direction_computation_first_order = std::make_shared<DirectionComputationProximalGradient>();
    direction_computation_first_order->addOptions(options, reporter);
    std::shared_ptr<DirectionComputation> direction_computation_second_order;
    direction_computation_second_order = std::make_shared<DirectionComputationProximalGradient>();
    direction_computation_second_order->addOptions(options, reporter);
    // ADD NEW DIRECTION COMPUTATION STRATEGIES HERE AND IN SWITCH BELOW //

    // Add options for line search strategies
    std::shared_ptr<LineSearch> line_search;
    line_search = std::make_shared<LineSearchBacktracking>();
    line_search->addOptions(options, reporter);
    // ADD NEW LINE SEARCH STRATEGIES HERE AND IN SWITCH BELOW //

    // Add options for parameter update strategies
    std::shared_ptr<ParameterUpdate> parameter_update;
    parameter_update = std::make_shared<ParameterUpdatePGStepsize>();
    parameter_update->addOptions(options, reporter);
    // ADD NEW PARAMETER UPDATE STRATEGIES HERE AND IN SWITCH BELOW //

}  // end addOptions

// Set options
void Strategies::getOptions(const Options* options, const Reporter* reporter)
{
    // Declare strategy names
    std::string space_partition_name;
    std::string direction_computation_first_order_name;
    std::string direction_computation_second_order_name;
    std::string line_search_name;

    // Read integer options
    options->valueAsString(reporter, "space_partition", space_partition_name);
    options->valueAsString(reporter, "line_search", line_search_name);
    options->valueAsString(reporter, "direction_computation_first_order",
                           direction_computation_first_order_name);
    options->valueAsString(reporter, "direction_computation_second_order",
                           direction_computation_first_order_name);
    options->valueAsString(reporter, "line_search", line_search_name);
    // Set space patiton strategy
    if (space_partition_name.compare("FirstOrderPartition") == 0)
    {
        space_partition_ = std::make_shared<SpacePartitionFirstOrder>();
    }
    else
    {
        space_partition_ = std::make_shared<SpacePartitionFirstOrder>();
    }

    // Set direction computation strategy
    if (direction_computation_first_order_name.compare("ProximalGradient") == 0)
    {
        direction_computation_first_order_ =
            std::make_shared<DirectionComputationProximalGradient>();
    }
    else
    {
        direction_computation_first_order_ =
            std::make_shared<DirectionComputationProximalGradient>();
    }

    if (direction_computation_first_order_name.compare("TruncatedNewton") == 0)
    {
        direction_computation_second_order_ =
            std::make_shared<DirectionComputationProximalGradient>();
    }
    else
    {
        direction_computation_second_order_ =
            std::make_shared<DirectionComputationProximalGradient>();
    }

    // Set line search strategy
    if (line_search_name.compare("Backtracking") == 0)
    {
        line_search_ = std::make_shared<LineSearchBacktracking>();
    }
    else
    {
        line_search_ = std::make_shared<LineSearchBacktracking>();
    }

    parameter_update_pg_stepsize_ = std::make_shared<ParameterUpdatePGStepsize>();

    // Set direction computation options
    space_partition_->getOptions(options, reporter);
    // Set direction computation options
    direction_computation_first_order_->getOptions(options, reporter);
    direction_computation_second_order_->getOptions(options, reporter);

    // Set line search options
    line_search_->getOptions(options, reporter);

    // set parameter update PG stepsize
    parameter_update_pg_stepsize_->getOptions(options, reporter);
}  // end getOptions

// Initialize
void Strategies::initialize(const Options* options, Quantities* quantities,
                            const Reporter* reporter)
{
    // Initialize space partition
    space_partition_->initialize(options, quantities, reporter);
    // Initialize direction computation
    direction_computation_first_order_->initialize(options, quantities, reporter);
    direction_computation_second_order_->initialize(options, quantities, reporter);

    // Initialize line search
    line_search_->initialize(options, quantities, reporter);
    // Initalize paramter update
    parameter_update_pg_stepsize_->initialize(options, quantities, reporter);
}  // end initialize

// Set iteration header
void Strategies::setIterationHeader()
{
    // Set iteration header string based on strategy objects
    iteration_header_ = "";
    if (direction_computation_first_order_->iterationHeader().length() > 0)
    {
        iteration_header_ += " ";
        iteration_header_ += direction_computation_first_order_->iterationHeader();
    }  // end if
    if (line_search_->iterationHeader().length() > 0)
    {
        iteration_header_ += " ";
        iteration_header_ += line_search_->iterationHeader();
    }  // end if

}  // end setIterationHeader

// Print header
void Strategies::printHeader(const Reporter* reporter)
{
    // Print header
    reporter->printf(
        R_SOLVER, R_BASIC,
        "\nSpace Partition strategy ...................... : %s\n"
        "Direction computation strategy (1st order)..... : %s\n"
        "Direction computation strategy (2nd order)..... : %s\n"
        "Line search strategy........................... : %s\n",
        space_partition_->name().c_str(), direction_computation_first_order_->name().c_str(),
        direction_computation_second_order_->name().c_str(), line_search_->name().c_str());

}  // end printHeader

// Print footer
void Strategies::printFooter(const Reporter* reporter) {}

}  // namespace FaRSA
