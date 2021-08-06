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
void DirectionComputationProximalGradient::addOptions(Options*        options,
                                                      const Reporter* reporter)
{
}  // end addOptions

// Set options
void DirectionComputationProximalGradient::getOptions(const Options*  options,
                                                      const Reporter* reporter)
{
}  // end getOptions

// Initialize
void DirectionComputationProximalGradient::initialize(const Options* options,
                                                      Quantities*    quantities,
                                                      const Reporter* reporter)
{
}

// Iteration header
std::string DirectionComputationProximalGradient::iterationHeader()
{
    return "  |Step| ";
}

// Iteration null values string
std::string DirectionComputationProximalGradient::iterationNullValues()
{
    return "---------";
}

// Compute direction
void DirectionComputationProximalGradient::computeDirection(
    const Options* options, Quantities* quantities, const Reporter* reporter,
    Strategies* strategies)
{
    // Initialize values
    setStatus(DC_UNSET);
    quantities->setTrialIterateToCurrentIterate();

    // try direction computation, terminate on any exception
    try
    {
        // Initialize boolean for evaluation
        bool evaluation_success = false;

        // Evaluate current objective
        evaluation_success =
            quantities->currentIterate()->evaluateObjective(*quantities);

        // Check for successful evaluation
        if (!evaluation_success)
        {
            THROW_EXCEPTION(DC_EVALUATION_FAILURE_EXCEPTION,
                            "Direction computation unsuccessful. "
                            "Objective Evaluation failed.");
        }

        // Evaluate current gradient
        evaluation_success =
            quantities->currentIterate()->evaluateGradient(*quantities);

        // Check for successful evaluation
        if (!evaluation_success)
        {
            THROW_EXCEPTION(
                DC_EVALUATION_FAILURE_EXCEPTION,
                "Direction computation unsuccessful. Evaluation failed.");
        }

        // Set trial iterate
        quantities->direction()->copy(
            *quantities->currentIterate()->gradient());
        quantities->direction()->scale(-1.0);

        // Set status
        setStatus(DC_SUCCESS);

        // Check for success
        THROW_EXCEPTION(DC_SUCCESS_EXCEPTION,
                        "Direction computation successful.")

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
    reporter->printf(R_SOLVER, R_PER_ITERATION, " %+.2e",
                     quantities->direction()->normInf());

}  // end computeDirection

}  // namespace FaRSA
