// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#include "FaRSAStrategies.hpp"

#include <stdlib.h>  // for terminateing the program

#include "FaRSADirectionComputationProximalGradient.hpp"
#include "FaRSADirectionComputationTruncatedNewton.hpp"
#include "FaRSALineSearchBacktracking.hpp"
#include "FaRSAParameterUpdate.hpp"
#include "FaRSAParameterUpdatePGStepsize.hpp"
#include "FaRSASpacePartitionFirstOrder.hpp"
#include "FaRSASpacePartitionGroupL1PGBased.hpp"

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
    options->addStringOption(reporter, "parameter_updates", "PGstepsize",
                             "A set of parameter update strategies to use. Sepearate by \";\"\n"
                             "Sample usages: PGstepsize;ProjectionL2BallRadius"
                             "Default     : PGstepsize.");
    // Add bool options
    options->addBoolOption(reporter, "strategies_verbose", true,
                           "A parameter controls whether should print more details.\n"
                           "Default     : true.");

    // Add options for space partition strategies
    std::shared_ptr<SpacePartition> space_partition;
    space_partition = std::make_shared<SpacePartitionFirstOrder>();
    space_partition->addOptions(options, reporter);
    space_partition = std::make_shared<SpacePartitionGroupL1PGBased>();
    space_partition->addOptions(options, reporter);
    // ADD NEW SPACE PARTITION STRATEGIES HERE AND IN SWITCH BELOW //

    // Add options for direction computation strategies
    std::shared_ptr<DirectionComputation> direction_computation_first_order;
    direction_computation_first_order = std::make_shared<DirectionComputationProximalGradient>();
    direction_computation_first_order->addOptions(options, reporter);
    std::shared_ptr<DirectionComputation> direction_computation_second_order;
    direction_computation_second_order = std::make_shared<DirectionComputationTruncatedNewton>();
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
    // get verbose_ values
    options->valueAsBool(reporter, "strategies_verbose", verbose_);

    // Declare strategy names
    std::string space_partition_name;
    std::string direction_computation_first_order_name;
    std::string direction_computation_second_order_name;
    std::string line_search_name;
    std::string parameter_updates_names;

    // Read integer options
    options->valueAsString(reporter, "space_partition", space_partition_name);
    options->valueAsString(reporter, "direction_computation_first_order", direction_computation_first_order_name);
    options->valueAsString(reporter, "direction_computation_second_order", direction_computation_second_order_name);
    options->valueAsString(reporter, "line_search", line_search_name);
    options->valueAsString(reporter, "parameter_updates", parameter_updates_names);

    // Set space patiton strategy
    if (space_partition_name.compare("FirstOrderPartition") == 0)
    {
        space_partition_ = std::make_shared<SpacePartitionFirstOrder>();
    }
    else if (space_partition_name.compare("GroupL1PGBased") == 0)
    {
        space_partition_ = std::make_shared<SpacePartitionGroupL1PGBased>();
    }
    else
    {
        reporter->printf(R_SOLVER, R_BASIC,
                         "Attempted to add un-recognized space partition strategy \"%s\". Terminate the program.\n",
                         space_partition_name.c_str());
        exit(EXIT_FAILURE);
    }
    // Set space partition computation options
    space_partition_->getOptions(options, reporter);

    // Set direction computation strategy
    if (direction_computation_first_order_name.compare("ProximalGradient") == 0)
    {
        direction_computation_first_order_ = std::make_shared<DirectionComputationProximalGradient>();
    }
    else
    {
        reporter->printf(
            R_SOLVER, R_BASIC,
            "Attempted to add un-recognized direction computation strategy \"%s\". Terminate the program.\n",
            direction_computation_first_order_name.c_str());
        exit(EXIT_FAILURE);
    }
    // Set direction computation options
    direction_computation_first_order_->getOptions(options, reporter);

    if (direction_computation_second_order_name.compare("TruncatedNewton") == 0)
    {
        direction_computation_second_order_ = std::make_shared<DirectionComputationTruncatedNewton>();
    }
    else
    {
        reporter->printf(
            R_SOLVER, R_BASIC,
            "Attempted to add un-recognized direction computation strategy \"%s\". Terminate the program.\n",
            direction_computation_second_order_name.c_str());
        exit(EXIT_FAILURE);
    }
    // Set direction computation options
    direction_computation_second_order_->getOptions(options, reporter);

    // Set line search strategy
    if (line_search_name.compare("Backtracking") == 0)
    {
        line_search_ = std::make_shared<LineSearchBacktracking>();
    }
    else
    {
        reporter->printf(R_SOLVER, R_BASIC,
                         "Attempted to add un-recognized line search strategy \"%s\". Terminate the program.\n",
                         line_search_name.c_str());
        exit(EXIT_FAILURE);
    }
    // Set line search options
    line_search_->getOptions(options, reporter);

    // Set a collection of parameter update strategies
    // parse the parameter_updates_set_of_names
    std::vector<std::string> parameter_update_lst;
    std::stringstream        s_stream(parameter_updates_names);
    while (s_stream.good())
    {
        std::string parameter_update;
        // get first string delimited by ;
        std::getline(s_stream, parameter_update, ';');
        parameter_update_lst.push_back(parameter_update);
    }
    parameter_updates_ = std::make_shared<ParameterUpdates>();
    for (std::string parameter_update : parameter_update_lst)
    {
        if (parameter_update.compare("PGstepsize") == 0)
        {
            auto temp = std::make_shared<ParameterUpdatePGStepsize>();
            // set parameter update PG stepsize
            temp->getOptions(options, reporter);
            bool is_success = parameter_updates_->add(temp);
            if (!is_success)
            {
                reporter->printf(R_SOLVER, R_BASIC,
                                 "Attempted to add parameter update strategy \"%s\", but failed. Terminate the "
                                 "program, which is not supposed to happen.\n",
                                 parameter_update.c_str());
                exit(EXIT_FAILURE);
            }
        }
        else
        {
            reporter->printf(
                R_SOLVER, R_BASIC,
                "Attempted to add un-recognized parameter update strategy \"%s\". Terminate the program.\n",
                parameter_update.c_str());
            exit(EXIT_FAILURE);
        }
    }

}  // end getOptions

// Initialize
void Strategies::initialize(const Options* options, Quantities* quantities, const Reporter* reporter)
{
    // Initialize space partition
    space_partition_->initialize(options, quantities, reporter);
    // Initialize direction computation
    direction_computation_first_order_->initialize(options, quantities, reporter);
    direction_computation_second_order_->initialize(options, quantities, reporter);

    // Initialize line search
    line_search_->initialize(options, quantities, reporter);
    // Initalize paramter updates
    parameter_updates_->initialize(options, quantities, reporter);
}  // end initialize

// Set iteration header
void Strategies::setIterationHeader()
{
    // Set iteration header string based on strategy objects
    iteration_header_ = "";
    if (space_partition_->iterationHeader().length() > 0)
    {
        iteration_header_ += " ";
        iteration_header_ += space_partition_->iterationHeader();
    }  // end if
    if (direction_computation_first_order_->iterationHeader().length() > 0)
    {
        iteration_header_ += " ";
        iteration_header_ += direction_computation_first_order_->iterationHeader();
    }  // end if
    if (direction_computation_second_order_->iterationHeader().length() > 0)
    {
        iteration_header_ += " ";
        iteration_header_ += direction_computation_second_order_->iterationHeader();
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
    bool use_second_order_direction = (space_partition_->name().compare("FirstOrderPartition") == 0);
    use_second_order_direction = !use_second_order_direction;
    reporter->printf(R_SOLVER, R_BASIC, "\n*************** Strategies Summary **************\n");
    reporter->printf(
        R_SOLVER, R_BASIC,
        "Space Partition strategy ...................... : %s\n"
        "Direction computation strategy (1st order)..... : %s\n"
        "Direction computation strategy (2nd order)..... : %s\n"
        "Line search strategy........................... : %s\n"
        "Parameter update strategies: .................. : %d strategy(ies) used.\n",
        space_partition_->name().c_str(), direction_computation_first_order_->name().c_str(),
        use_second_order_direction ? direction_computation_second_order_->name().c_str() : "Not applicable",
        line_search_->name().c_str(), parameter_updates_->size());
    if (verbose_)
    {
        reporter->printf(R_SOLVER, R_BASIC,
                         "\n------------------ Deatils ----------------\n"
                         "Line search strategy:\n%s"
                         "Parameter update strategies:\n%s",
                         line_search_->details().c_str(), parameter_updates_->details().c_str());
    }

}  // end printHeader

// Print footer
void Strategies::printFooter(const Reporter* reporter) {}

}  // namespace FaRSA
