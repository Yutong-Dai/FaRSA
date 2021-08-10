// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#include "FaRSASolver.hpp"

#include <cmath>
#include <iostream>

#include "FaRSADeclarations.hpp"
#include "FaRSADefinitions.hpp"
#include "FaRSAException.hpp"
#include "FaRSAVersion.hpp"

namespace FaRSA
{
// Constructor
FaRSASolver::FaRSASolver()
{
    // Declare stream report
    std::shared_ptr<StreamReport> s(new StreamReport("default", R_SOLVER, R_PER_ITERATION));

    // Set stream report to standard output
    s->setStream(&std::cout);

    // Add stream report to reporter
    reporter_.addReport(s);

    // Add options
    addOptions();

}  // end constructor

// Destructor
FaRSASolver::~FaRSASolver()
{
    // Delete reporter
    reporter_.deleteReports();

}  // end destructor

// Add options
void FaRSASolver::addOptions()
{
    // Add options for quantities
    quantities_.addOptions(&options_, &reporter_);

    // Add options for strategies
    strategies_.addOptions(&options_, &reporter_);
}  // end addOptions

// Set options
void FaRSASolver::getOptions()
{
    // Set quantities options
    quantities_.getOptions(&options_, &reporter_);

    // Set strategies options
    strategies_.getOptions(&options_, &reporter_);

}  // end getOptions

// Solution
void FaRSASolver::solution(double vector[])
{
    // Copy values to solution
    for (int i = 0; i < (int)quantities_.numberOfVariables(); i++)
    {
        vector[i] = quantities_.currentIterate()->vector()->values()[i];
    }

}  // end solution

// Optimize
void FaRSASolver::optimize(const std::shared_ptr<FunctionSmooth>    function_smooth,
                           const std::shared_ptr<FunctionNonsmooth> function_nonsmooth,
                           const Vector&                            initail_point)
{
    // Initialize solver status
    setStatus(FARSA_UNSET);

    // (Re)set options
    getOptions();

    // try to run algorithm, terminate on any exception
    try
    {
        // (Re)initialize quantities
        bool initialization_success =
            quantities_.initialize(function_smooth, function_nonsmooth, initail_point);

        // Check for initialization success
        if (!initialization_success)
        {
            THROW_EXCEPTION(FARSA_INITIALIZATION_FAILURE_EXCEPTION, "Initialization failed.");
        }
        /**
         * The scale is determined within the quantities iniatilization and the objective value
         * and the gradient are properly scaled, which are evaluated at the initial point
         */
        // // Evaluate all functions at current iterate
        // evaluateFunctionsAtCurrentIterate();

        // // Determine problem scaling
        // quantities_.currentIterate()->determineScale(quantities_);

        // // Scale evaluated objective
        // quantities_.currentIterate()->scaleObjective();

        // // Scale evaluated gradient
        // quantities_.currentIterate()->scaleGradient();

        // Store norm of initial point (for termination check)
        double initial_iterate_norm = quantities_.currentIterate()->vector()->norm2();
        // Store norm of proximal gradient step at the initial point (for termination check)
        computeProximalGradientUpdateAtCurrentIterate();
        double initial_proximal_gradient_step_norm =
            quantities_.currentIterate()->proximalGraidentStep()->normInf();
        // Initialize strategies
        strategies_.initialize(&options_, &quantities_, &reporter_);

        // Print header
        printHeader();

        // Set iteration header
        strategies_.setIterationHeader();

        // (Outer) Loop
        while (true)
        {
            // Print iteration header
            printIterationHeader();

            // Print quantities iteration values
            quantities_.printIterationValues(&reporter_);

            // Flush buffer
            reporter_.flushBuffer();

            // Check termination conditions
            computeProximalGradientUpdateAtCurrentIterate();

            if (quantities_.currentIterate()->proximalGraidentStep()->normInf() <=
                quantities_.stationarityTolerance() *
                    fmax(1.0, initial_proximal_gradient_step_norm))
            {
                THROW_EXCEPTION(FARSA_SUCCESS_EXCEPTION, "Stationary point found.");
            }
            if (quantities_.iterationCounter() >= quantities_.iterationLimit())
            {
                THROW_EXCEPTION(FARSA_ITERATION_LIMIT_EXCEPTION,
                                "Iteration limit has been reached.");
            }
            if ((clock() - quantities_.startTime()) / (double)CLOCKS_PER_SEC >=
                quantities_.cpuTimeLimit())
            {
                THROW_EXCEPTION(FARSA_CPU_TIME_LIMIT_EXCEPTION, "CPU time limit has been reached.");
            }
            if (quantities_.currentIterate()->vector()->norm2() >=
                quantities_.iterateNormTolerance() * fmax(1.0, initial_iterate_norm))
            {
                THROW_EXCEPTION(FARSA_ITERATE_NORM_LIMIT_EXCEPTION,
                                "Iterates appear to be diverging.");
            }
            // Perform space partition
            strategies_.spacePartition()->partitionSpace(&options_, &quantities_, &reporter_,
                                                         &strategies_);
            // Check status
            if (strategies_.spacePartition()->status() != SP_SUCCESS)
            {
                THROW_EXCEPTION(FARSA_SPACE_PARTITION_FAILURE_EXCEPTION, "Space partition failed.")
            }
            // Compute direction for the first order group of variables
            strategies_.directionComputationFirstOrder()->computeDirection(
                &options_, &quantities_, &reporter_, &strategies_);
            // Check status
            if (strategies_.directionComputationFirstOrder()->status() != DC_SUCCESS)
            {
                THROW_EXCEPTION(FARSA_DIRECTION_COMPUTATION_FAILURE_EXCEPTION,
                                "Direction computation for the first order group failed.");
            }
            // Compute direction for the second order group of variables
            strategies_.directionComputationSecondOrder()->computeDirection(
                &options_, &quantities_, &reporter_, &strategies_);
            // Check status
            if (strategies_.directionComputationFirstOrder()->status() != DC_SUCCESS)
            {
                THROW_EXCEPTION(FARSA_DIRECTION_COMPUTATION_FAILURE_EXCEPTION,
                                "Direction computation for the second order group failed.");
            }
            // quantities_.direction()->print(&reporter_, "\ndirection");
            // Run line search
            strategies_.lineSearch()->runLineSearch(&options_, &quantities_, &reporter_,
                                                    &strategies_);

            // Check status
            if (strategies_.lineSearch()->status() != LS_SUCCESS)
            {
                THROW_EXCEPTION(FARSA_LINE_SEARCH_FAILURE_EXCEPTION, "Line search failed.");
            }

            // Update Parameters (this has to be done before performing the iterate update)
            strategies_.parameterUpdatePGStepsize()->update(&options_, &quantities_, &reporter_,
                                                            &strategies_);
            // Check status
            if (strategies_.parameterUpdatePGStepsize()->status() != PU_SUCCESS)
            {
                THROW_EXCEPTION(FARSA_PARAMETER_UPDATE_FAILURE_EXCEPTION,
                                "Parameter Update for the proximal gradient computation failed.")
            }
            // Update iterate
            quantities_.setCurrentIterate(quantities_.trialIterate());

            // Increment iteration counter
            quantities_.incrementIterationCounter();

            // Evaluate all functions at current iterate
            evaluateFunctionstAtCurrentIterate();

            // Print end of line
            reporter_.printf(R_SOLVER, R_PER_ITERATION, "\n");

        }  // end while

    }  // end try

    // catch exceptions
    catch (FARSA_INITIALIZATION_FAILURE_EXCEPTION& exec)
    {
        setStatus(FARSA_INITIALIZATION_FAILURE);
    }
    catch (FARSA_SUCCESS_EXCEPTION& exec)
    {
        setStatus(FARSA_SUCCESS);
    }
    // limits exception
    catch (FARSA_ITERATION_LIMIT_EXCEPTION& exec)
    {
        setStatus(FARSA_ITERATION_LIMIT);
    }
    catch (FARSA_CPU_TIME_LIMIT_EXCEPTION& exec)
    {
        setStatus(FARSA_CPU_TIME_LIMIT);
    }
    catch (FARSA_ITERATE_NORM_LIMIT_EXCEPTION& exec)
    {
        setStatus(FARSA_ITERATE_NORM_LIMIT);
    }
    catch (FARSA_FUNCTION_EVALUATION_LIMIT_EXCEPTION& exec)
    {
        setStatus(FARSA_FUNCTION_EVALUATION_LIMIT);
    }
    catch (FARSA_GRADIENT_EVALUATION_LIMIT_EXCEPTION& exec)
    {
        setStatus(FARSA_GRADIENT_EVALUATION_LIMIT);
    }
    catch (FARSA_HESSIAN_VECTOR_PRODUCT_EVALUATION_LIMIT_EXCEPTION& exec)
    {
        setStatus(FARSA_HESSIAN_VECTOR_PRODUCT_EVALUATION_LIMIT);
    }

    // evaluation/computation failure
    catch (FARSA_FUNCTION_EVALUATION_FAILURE_EXCEPTION& exec)
    {
        setStatus(FARSA_FUNCTION_EVALUATION_FAILURE);
    }
    catch (FARSA_GRADIENT_EVALUATION_FAILURE_EXCEPTION& exec)
    {
        setStatus(FARSA_GRADIENT_EVALUATION_FAILURE);
    }
    catch (FARSA_PROXIMAL_GRADIENT_COMPUTATION_FAILURE_EXCEPTION& exec)
    {
        setStatus(FARSA_PROXIMAL_GRADIENT_COMPUTATION_FAILURE);
    }
    catch (FARSA_HESSIAN_VECTOR_PRODUCT_EVALUATION_FAILURE_EXCEPTION& exec)
    {
        setStatus(FARSA_HESSIAN_VECTOR_PRODUCT_EVALUATION_FAILURE);
    }
    catch (FARSA_FUNCTION_EVALUATION_ASSERT_EXCEPTION& exec)
    {
        setStatus(FARSA_FUNCTION_EVALUATION_ASSERT);
    }
    catch (FARSA_GRADIENT_EVALUATION_ASSERT_EXCEPTION& exec)
    {
        setStatus(FARSA_GRADIENT_EVALUATION_ASSERT);
    }
    // strategies exception
    catch (FARSA_SPACE_PARTITION_FAILURE_EXCEPTION& exec)
    {
        setStatus(FARSA_SPACE_PARTITION_FAILURE);
    }
    catch (FARSA_DIRECTION_COMPUTATION_FAILURE_EXCEPTION& exec)
    {
        setStatus(FARSA_DIRECTION_COMPUTATION_FAILURE);
    }
    catch (FARSA_LINE_SEARCH_FAILURE_EXCEPTION& exec)
    {
        setStatus(FARSA_LINE_SEARCH_FAILURE);
    }
    catch (FARSA_PARAMETER_UPDATE_FAILURE_EXCEPTION& exec)
    {
        setStatus(FARSA_PARAMETER_UPDATE_FAILURE);
    }

    // lowlevel-objects exception
    catch (FARSA_MATRIX_EXCEPTION& exec)
    {
        setStatus(FARSA_MATRIX);
    }
    catch (FARSA_MATRIX_ASSERT_EXCEPTION& exec)
    {
        setStatus(FARSA_MATRIX_ASSERT);
    }
    catch (FARSA_VECTOR_EXCEPTION& exec)
    {
        setStatus(FARSA_VECTOR);
    }
    catch (FARSA_VECTOR_ASSERT_EXCEPTION& exec)
    {
        setStatus(FARSA_VECTOR_ASSERT);
    }

    // Print end of line
    reporter_.printf(R_SOLVER, R_PER_ITERATION, "\n");

    // // Check whether to finalize problem solution
    // if (status() != FARSA_INITIALIZATION_FAILURE && status() != FARSA_FUNCTION_EVALUATION_FAILURE
    // &&
    //     status() != FARSA_GRADIENT_EVALUATION_FAILURE)
    // {
    //     // Finalize problem solution
    //     problem->finalizeSolution(quantities_.currentIterate()->vector()->values(),
    //                               quantities_.currentIterate()->objective(),
    //                               quantities_.currentIterate()->gradient()->values());
    // }

    // Finalize
    quantities_.finalize();

    // Print footer
    printFooter();

}  // end optimize

// Evaluate smooth and nonsmooth functions at current iterate
void FaRSASolver::evaluateFunctionstAtCurrentIterate()
{
    // Evaluate objective
    bool evaluation_success = quantities_.currentIterate()->evaluateObjectiveAll(quantities_);

    // Check for evaluation success
    if (!evaluation_success)
    {
        THROW_EXCEPTION(FARSA_FUNCTION_EVALUATION_FAILURE_EXCEPTION, "Function evaluation failed.");
    }

}  // end evaluateFunctionsAtCurrentIterate

void FaRSASolver::computeProximalGradientUpdateAtCurrentIterate()
{
    bool evaluation_success = quantities_.currentIterate()->evaluateObjectiveSmooth(quantities_);
    // Check for evaluation success
    if (!evaluation_success)
    {
        THROW_EXCEPTION(FARSA_FUNCTION_EVALUATION_FAILURE_EXCEPTION,
                        "Function smooth evaluation failed.");
    }
    evaluation_success = quantities_.currentIterate()->evaluateGradientSmooth(quantities_);
    // Check for evaluation success
    if (!evaluation_success)
    {
        THROW_EXCEPTION(FARSA_FUNCTION_EVALUATION_FAILURE_EXCEPTION,
                        "Graient smooth evaluation failed.");
    }
    bool computation_success =
        quantities_.currentIterate()->computeProximalGradientUpdate(quantities_);
    // Check for computation success
    if (!computation_success)
    {
        THROW_EXCEPTION(FARSA_PROXIMAL_GRADIENT_COMPUTATION_FAILURE_EXCEPTION,
                        "Proximal gradient computation failed.");
    }
}

// Print footer
void FaRSASolver::printFooter()
{
    // Print footer
    reporter_.printf(R_SOLVER, R_BASIC, "\nEXIT: ");

    // Print exit status
    switch (status())
    {
        case FARSA_UNSET:
            reporter_.printf(R_SOLVER, R_BASIC,
                             "Exit status wasn't set! This wasn't supposed to happen!");
            break;
        case FARSA_SUCCESS:
            reporter_.printf(R_SOLVER, R_BASIC, "Stationary point found.");
            break;
        case FARSA_CPU_TIME_LIMIT:
            reporter_.printf(R_SOLVER, R_BASIC, "CPU time limit reached.");
            break;
        case FARSA_ITERATE_NORM_LIMIT:
            reporter_.printf(R_SOLVER, R_BASIC, "Iterates seem to be diverging.");
            break;
        case FARSA_ITERATION_LIMIT:
            reporter_.printf(R_SOLVER, R_BASIC, "Iteration limit reached.");
            break;
        case FARSA_FUNCTION_EVALUATION_LIMIT:
            reporter_.printf(R_SOLVER, R_BASIC, "Function evaluation limit reached.");
            break;
        case FARSA_GRADIENT_EVALUATION_LIMIT:
            reporter_.printf(R_SOLVER, R_BASIC, "Gradient evaluation limit reached.");
            break;
        case FARSA_HESSIAN_VECTOR_PRODUCT_EVALUATION_LIMIT:
            reporter_.printf(R_SOLVER, R_BASIC, "Hessian-vector product evaluation limit reached.");
            break;
        case FARSA_INITIALIZATION_FAILURE:
            reporter_.printf(R_SOLVER, R_BASIC,
                             "Initialization failure! Check definition of problem.");
            break;
        case FARSA_FUNCTION_EVALUATION_FAILURE:
            reporter_.printf(R_SOLVER, R_BASIC,
                             "Function evaluation failure! Check definition of function.");
            break;
        case FARSA_GRADIENT_EVALUATION_FAILURE:
            reporter_.printf(R_SOLVER, R_BASIC,
                             "Gradient evaluation failure! Check definition of function.");
            break;
        case FARSA_PROXIMAL_GRADIENT_COMPUTATION_FAILURE:
            reporter_.printf(R_SOLVER, R_BASIC,
                             "Proximal graident computation failure! Check definition of nonsmooth "
                             "function proximalGradientUpdate method.");
            break;
        case FARSA_HESSIAN_VECTOR_PRODUCT_EVALUATION_FAILURE:
            reporter_.printf(R_SOLVER, R_BASIC,
                             "Hessian-vector product failure! Check definition of function.");
            break;
        case FARSA_FUNCTION_EVALUATION_ASSERT:
            reporter_.printf(R_SOLVER, R_BASIC,
                             "Function evaluation assert failure! This wasn't supposed to happen!");
            break;
        case FARSA_GRADIENT_EVALUATION_ASSERT:
            reporter_.printf(R_SOLVER, R_BASIC,
                             "Gradient evaluation assert failure! This wasn't supposed to happen!");
            break;
        case FARSA_SPACE_PARTITION_FAILURE:
            reporter_.printf(R_SOLVER, R_BASIC, "Space Decomposition failure.");
            break;
        case FARSA_DIRECTION_COMPUTATION_FAILURE:
            reporter_.printf(R_SOLVER, R_BASIC, "Direction computation failure.");
            break;
        case FARSA_LINE_SEARCH_FAILURE:
            reporter_.printf(R_SOLVER, R_BASIC, "Line search failure.");
            break;
        case FARSA_PARAMETER_UPDATE_FAILURE:
            reporter_.printf(R_SOLVER, R_BASIC, "Parameter update failure.");
            break;
        default:
            reporter_.printf(R_SOLVER, R_BASIC,
                             "Unknown exit status! This wasn't supposed to happen!");
            break;
    }  // end switch

    // Check whether to print final data
    if (status() != FARSA_INITIALIZATION_FAILURE && status() != FARSA_FUNCTION_EVALUATION_FAILURE)
    {
        // Print quantities footer
        quantities_.printFooter(&reporter_);

        // Print strategies footer
        strategies_.printFooter(&reporter_);

    }  // end if

}  // end printFooter

// Print header
void FaRSASolver::printHeader()
{
    // Print header
    reporter_.printf(R_SOLVER, R_BASIC,
                     "+-------------------------------------------------------------+\n"
                     "|             FaRSA = Fast Reduced Space Algorithm            |\n"
                     "| FaRSA is released as open source code under the ??? License |\n"
                     "+-------------------------------------------------------------+\n"
                     "\n"
                     "This is FaRSA version %s\n"
                     "\n",
                     FARSA_VERSION);

    // Print quantities header
    quantities_.printHeader(&reporter_);

    // Print strategies header
    strategies_.printHeader(&reporter_);

}  // end printHeader

// Print iteration header
void FaRSASolver::printIterationHeader()
{
    if (quantities_.iterationCounter() == 0)
    {
        reporter_.printf(R_SOLVER, R_PER_ITERATION, "\n");
    }
    if (quantities_.iterationCounter() % 20 == 0)
    {
        std::string b(
            quantities_.iterationHeader().length() + strategies_.iterationHeader().length(), '-');
        reporter_.printf(
            R_SOLVER, R_PER_ITERATION, "%s\n",
            (b + "\n" + quantities_.iterationHeader() + strategies_.iterationHeader() + "\n" + b)
                .c_str());
    }  // end if

}  // end printIterationHeader

}  // namespace FaRSA
