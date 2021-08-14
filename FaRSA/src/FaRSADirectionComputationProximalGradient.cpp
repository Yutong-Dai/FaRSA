// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#include "FaRSADirectionComputationProximalGradient.hpp"

#include "FaRSADeclarations.hpp"
#include "FaRSADefinitions.hpp"

namespace FaRSA
{
// Add options
void DirectionComputationProximalGradient::addOptions(Options* options, const Reporter* reporter) {}  // end addOptions

// Set options
void DirectionComputationProximalGradient::getOptions(const Options* options, const Reporter* reporter) {
}  // end getOptions

// Initialize
void DirectionComputationProximalGradient::initialize(const Options* options, Quantities* quantities,
                                                      const Reporter* reporter)
{
}

// Iteration header
std::string DirectionComputationProximalGradient::iterationHeader() { return "|Direction|"; }

// Iteration null values string
std::string DirectionComputationProximalGradient::iterationNullValues() { return "---------"; }

// Compute direction
void DirectionComputationProximalGradient::computeDirection(const Options* options, Quantities* quantities,
                                                            const Reporter* reporter, Strategies* strategies)
{
    // Initialize values
    setStatus(DC_UNSET);
    // if perform computation
    if (performCompuation())
    {
        // try direction computation, terminate on any exception
        try
        {
            // Initialize boolean for evaluation
            bool evaluation_success = false;

            // Evaluate current objective
            evaluation_success = quantities->currentIterate()->evaluateObjectiveAll(*quantities);

            // Check for successful evaluation
            if (!evaluation_success)
            {
                THROW_EXCEPTION(DC_EVALUATION_FAILURE_EXCEPTION,
                                "Direction computation unsuccessful. "
                                "Objective Smooth Evaluation failed.");
            }
            // Evaluate current gradient
            evaluation_success = quantities->currentIterate()->evaluateGradientSmooth(*quantities);

            // Check for successful evaluation
            if (!evaluation_success)
            {
                THROW_EXCEPTION(DC_EVALUATION_FAILURE_EXCEPTION,
                                "Direction computation unsuccessful. Gradient "
                                "Smooth Evaluation failed.");
            }
            // Compute proximal graident
            evaluation_success = quantities->currentIterate()->computeProximalGradientUpdate(*quantities);
            if (!evaluation_success)
            {
                THROW_EXCEPTION(DC_EVALUATION_FAILURE_EXCEPTION,
                                "Direction computation unsuccessful. Proximal "
                                "gradient computation failed.");
            }
            // Set direction possibly in a low dimension;  then 0 is paading to
            // non-working indicies
            // directon = proximalGradientUpdate - currentIterate
            Vector search_direction_actual(quantities->numberOfVariables());
            auto   step_ptr = quantities->currentIterate()->proximalGraidentStep();
            for (auto i : *(quantities->indiciesWorking()))
            {
                search_direction_actual.valuesModifiable()[i] = (*step_ptr).values()[i];
            }
            quantities->directionFirstOrder()->copy(search_direction_actual);

            // Set status
            setStatus(DC_SUCCESS);
            // Check for success
            THROW_EXCEPTION(DC_SUCCESS_EXCEPTION, "Direction computation successful.")

        }  // end try

        // catch exceptions
        catch (DC_SUCCESS_EXCEPTION& exec)
        {
            setStatus(DC_SUCCESS);
        }
        catch (DC_EVALUATION_FAILURE_EXCEPTION& exec)
        {
            setStatus(DC_EVALUATION_FAILURE);
        }
        // Print iteration information
        // reporter->printf(R_SOLVER, R_PER_ITERATION, " %+.2e", quantities->direction()->normInf());
    }
    else
    {
        setStatus(DC_SKIPPED);
    }

}  // end computeDirection

}  // namespace FaRSA
