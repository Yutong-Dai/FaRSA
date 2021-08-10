// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#include "FaRSALineSearchBacktracking.hpp"

#include <cmath>
#include <iostream>

#include "FaRSADeclarations.hpp"
#include "FaRSADefinitions.hpp"

namespace FaRSA
{
// Add options
void LineSearchBacktracking::addOptions(Options* options, const Reporter* reporter)
{
    // Add bool options
    options->addBoolOption(reporter, "LSB_fail_on_small_stepsize", false,
                           "Indicator for whether to indicate failure on small stepsize.\n"
                           "Default     : false.");

    // Add double options
    options->addDoubleOption(reporter, "LSB_stepsize_initial", 1.0, 0.0, FARSA_DOUBLE_INFINITY,
                             "Initial stepsize to be used in the first iteration.  Note that\n"
                             "              the initial stepsize used in the line search in "
                             "subsequent\n"
                             "              iterations is set the minimum of this value and a "
                             "factor times\n"
                             "              the stepsize accepted in the previous iteration.\n"
                             "Default     : 1.0.");
    options->addDoubleOption(reporter, "LSB_stepsize_minimum", 1e-20, 0.0, FARSA_DOUBLE_INFINITY,
                             "Tolerance for determining an insufficient stepsize.  If the\n"
                             "              line search yields a stepsize below this tolerance, "
                             "then the\n"
                             "              algorithm may terminate with a message of a small "
                             "stepsize.\n"
                             "Default     : 1e-20.");
    options->addDoubleOption(reporter, "LSB_stepsize_sufficient_decrease_threshold", 1e-10, 0.0,
                             1.0,
                             "Sufficient decrease constant for the Armijo line search.\n"
                             "Default     : 1e-10.");
    options->addDoubleOption(reporter, "LSB_stepsize_sufficient_decrease_fudge_factor", 1e-10, 0.0,
                             FARSA_DOUBLE_INFINITY,
                             "Sufficient decrease fudge factor.\n"
                             "Default     : 1e-10.");
    options->addDoubleOption(reporter, "LSB_stepsize_decrease_factor", 5e-01, 0.0, 1.0,
                             "Factor for updating the stepsize during the line search.\n"
                             "Default     : 5e-01.");
    // Add integer options
    // options->addIntegerOption(reporter, "linesearch_max_backtrack", 100, 0, FARSA_INT_INFINITY,
    //                           "Limit on the number of backtrack-linesearch will be performed.\n"
    //                           "Default     : 100.");
}  // end addOptions

// Set options
void LineSearchBacktracking::getOptions(const Options* options, const Reporter* reporter)
{
    // Read bool options
    options->valueAsBool(reporter, "LSB_fail_on_small_stepsize", fail_on_small_stepsize_);

    // Read options
    options->valueAsDouble(reporter, "LSB_stepsize_initial", stepsize_initial_);
    options->valueAsDouble(reporter, "LSB_stepsize_minimum", stepsize_minimum_);
    options->valueAsDouble(reporter, "LSB_stepsize_sufficient_decrease_threshold",
                           stepsize_sufficient_decrease_threshold_);
    options->valueAsDouble(reporter, "LSB_stepsize_sufficient_decrease_fudge_factor",
                           stepsize_sufficient_decrease_fudge_factor_);
    options->valueAsDouble(reporter, "LSB_stepsize_decrease_factor", stepsize_decrease_factor_);
    // options->valueAsInteger(reporter, "linesearch_max_backtrack", max_backtrack_);

}  // end getOptions

// Initialize
void LineSearchBacktracking::initialize(const Options* options, Quantities* quantities,
                                        const Reporter* reporter)
{
    quantities->setStepsizeLineSearch(fmax(stepsize_minimum_, stepsize_initial_));
}

// Run line search - get the acceptable trial iterate; move from current iterate to the trial
// iterate is done in the solver
void LineSearchBacktracking::runLineSearch(const Options* options, Quantities* quantities,
                                           const Reporter* reporter, Strategies* strategies)
{
    // Initialize values
    setStatus(LS_UNSET);
    // set number_of_backtrack_ to zero;
    reset();
    // try line search, terminate on any exception
    try
    {
        // Evaluate composite objective at current point
        bool evaluation_success = quantities->currentIterate()->evaluateObjectiveAll(*quantities);

        // Check for successful evaluation
        if (!evaluation_success)
        {
            quantities->setStepsizeLineSearch(0.0);
            THROW_EXCEPTION(LS_EVALUATION_FAILURE_EXCEPTION,
                            "Line search unsuccessful. Evaluation of Objective failed.");
        }

        // Initialize stepsize
        quantities->setStepsizeLineSearch(fmax(stepsize_minimum_, stepsize_initial_));

        // line search along the proximal gradient direction
        double directional_derivative = 0.0;
        if (quantities->directionType() == DC_PROXIMAL_GRADIENT)
        {
            // Compute directional derivative upper-bound
            auto   direction = quantities->direction();
            double directional_derivative = direction->innerProduct(*direction);
            directional_derivative /= -quantities->stepsizeProximalGradient();
        }
        else
        {
            THROW_EXCEPTION(LS_EVALUATION_FAILURE_EXCEPTION, "Unsupported Direction Type.")
        }

        // Loop
        while (true)
        {
            // Declare new point
            quantities->setTrialIterate(quantities->currentIterate()->makeNewLinearCombination(
                1.0, quantities->stepsizeLineSearch(), *quantities->direction()));

            // Evaluate trial objective for both smooth and nonsmooth function
            evaluation_success = quantities->trialIterate()->evaluateObjectiveAll(*quantities);

            // Check for successful evaluation
            if (evaluation_success)
            {
                // Check for sufficient decrease
                bool sufficient_decrease =
                    (quantities->trialIterate()->objectiveAll() -
                         quantities->currentIterate()->objectiveAll() <=
                     stepsize_sufficient_decrease_threshold_ * quantities->stepsizeLineSearch() *
                         directional_derivative);
                // auto LHS = quantities->trialIterate()->objectiveAll() -
                //            quantities->currentIterate()->objectiveAll();
                // auto RHS = stepsize_sufficient_decrease_threshold_ *
                //            quantities->stepsizeLineSearch() * directional_derivative;
                // std::cout << "\n bak:" << number_of_backtrack_ << " LHS : " << LHS
                //           << " RHS : " << RHS << " LHS<=RHS: " << sufficient_decrease <<
                //           std::endl;
                // Check Armijo condition
                if (sufficient_decrease)
                {
                    THROW_EXCEPTION(LS_SUCCESS_EXCEPTION, "Line search successful.");
                }  // end if

            }  // end if

            // Check if stepsize below minimum
            if (quantities->stepsizeLineSearch() <= stepsize_minimum_)
            {
                // Check for failure on small stepsize
                if (fail_on_small_stepsize_)
                {
                    THROW_EXCEPTION(LS_STEPSIZE_TOO_SMALL_EXCEPTION,
                                    "Line search unsuccessful.  Stepsize too small.");
                }

                // Evaluate objective at trial iterate
                evaluation_success = quantities->trialIterate()->evaluateObjectiveAll(*quantities);

                // Check for evaluation success
                if (evaluation_success)
                {
                    // Check for decrease
                    if (quantities->trialIterate()->objectiveAll() <
                        quantities->currentIterate()->objectiveAll())
                    {
                        THROW_EXCEPTION(LS_SUCCESS_EXCEPTION, "Line search successful.");

                    }  // end if

                }  // end if

                // Set null step
                quantities->setStepsizeLineSearch(0.0);

                // Set new point
                quantities->setTrialIterateToCurrentIterate();

                // Terminate
                THROW_EXCEPTION(LS_SUCCESS_EXCEPTION, "Line search successful.");

            }  // end if

            // Update stepsize
            quantities->setStepsizeLineSearch(stepsize_decrease_factor_ *
                                              quantities->stepsizeLineSearch());
            // begin to backtrack
            number_of_backtrack_ += 1;

            // if (number_of_backtrack_ > max_backtrack_)
            // {
            //     THROW_EXCEPTION(LS_MAXBACKTRACK_LIMIT_EXCEPTION,
            //                     "Line search unsuccessful.  Reach the maximum backtrack limit.");
            // }

        }  // end while

    }  // end try

    // catch exceptions
    catch (LS_SUCCESS_EXCEPTION& exec)
    {
        setStatus(LS_SUCCESS);
    }
    catch (LS_EVALUATION_FAILURE_EXCEPTION& exec)
    {
        setStatus(LS_EVALUATION_FAILURE);
    }
    catch (LS_STEPSIZE_TOO_SMALL_EXCEPTION& exec)
    {
        setStatus(LS_STEPSIZE_TOO_SMALL);
    }  // end catch

    // Print iteration information
    reporter->printf(R_SOLVER, R_PER_ITERATION, " %+.2e", quantities->stepsizeLineSearch());

}  // end runLineSearch

}  // namespace FaRSA
