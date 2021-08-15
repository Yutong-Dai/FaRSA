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
#include "FaRSASpacePartitionGroupL1PGBased.hpp"

namespace FaRSA
{
// Add options
void LineSearchBacktracking::addOptions(Options* options, const Reporter* reporter)
{
    // Add bool options
    options->addBoolOption(reporter, "LSB_fail_on_small_stepsize", false,
                           "Indicator for whether to indicate failure on small stepsize.\n"
                           "Default     : false.");
    options->addBoolOption(reporter, "LSB_verbose", true,
                           "A parameter controls whether should print more details.\n"
                           "Default     : true.");

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
    options->addDoubleOption(reporter, "LSB_stepsize_sufficient_decrease_threshold", 1e-10, 0.0, 1.0,
                             "Sufficient decrease constant for the Armijo line search.\n"
                             "Default     : 1e-10.");
    options->addDoubleOption(reporter, "LSB_stepsize_sufficient_decrease_fudge_factor", 1e-10, 0.0,
                             FARSA_DOUBLE_INFINITY,
                             "Sufficient decrease fudge factor.\n"
                             "Default     : 1e-10.");
    options->addDoubleOption(reporter, "LSB_stepsize_decrease_factor", 5e-01, 0.0, 1.0,
                             "Factor for updating the stepsize during the line search.\n"
                             "Default     : 5e-01.");
    // Add string options
    options->addStringOption(
        reporter, "LSB_direction_search_type", "pure",
        "Type used to construct the search direction from the first order search direction and the second order search "
        "direction. Currently supported values: \"pure\"/\"hybrid\".\n"
        "If it is set to \"pure\", then the search direction is constructed by using only either the first order "
        "search direction or the second order search direction.\n"
        "If it is set to \"hybrid\", then the search direction is constructed by using the both directions.\n"
        "Default     : pure");
    options->addStringOption(reporter, "LSB_directional_derivative_first_order_type", "proximal_step",
                             "(Upper bound of the) Directional derivative along the first order search direction.\n"
                             "Currently supported values: proximal_step.\n"
                             "Default     : proximal_step");
    options->addStringOption(reporter, "LSB_directional_derivative_second_order_type", "gradient",
                             "(Upper bound of) Directional derivative along the second order search direction.\n"
                             "Currently supported values: gradient.\n"
                             "Default     : gradient");
    options->addStringOption(reporter, "LSB_method_pure_first_order", "armijo",
                             "Method used to perform backtrack linesearch along the first order search direction. This "
                             "requires \"LSB_direction_search_type\" to be set to \"pure\"\n"
                             "Currently supported values: armijo.\n"
                             "Default     : armijo");
    options->addStringOption(
        reporter, "LSB_method_pure_second_order", "projected_armijo_groupl1",
        "Method used to perform backtrack linesearch along the second order search direction. This "
        "requires \"LSB_direction_search_type\" to be set to \"pure\"\n"
        "Currently supported values: projected_armijo_groupl1.\n"
        "Default     : projected_armijo_groupl1");
    options->addStringOption(reporter, "LSB_method_hybrid", "armijo",
                             "Method used to perform backtrack linesearch along the hybrid search direction. This "
                             "requires \"LSB_direction_search_type\" to be set to \"hybrid\"\n"
                             "Currently supported values: dogleg.\n"
                             "Default     : dogleg");
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
    options->valueAsBool(reporter, "LSB_verbose", verbose_);
    // Read options
    options->valueAsDouble(reporter, "LSB_stepsize_initial", stepsize_initial_);
    options->valueAsDouble(reporter, "LSB_stepsize_minimum", stepsize_minimum_);
    options->valueAsDouble(reporter, "LSB_stepsize_sufficient_decrease_threshold",
                           stepsize_sufficient_decrease_threshold_);
    options->valueAsDouble(reporter, "LSB_stepsize_sufficient_decrease_fudge_factor",
                           stepsize_sufficient_decrease_fudge_factor_);
    options->valueAsDouble(reporter, "LSB_stepsize_decrease_factor", stepsize_decrease_factor_);
    options->valueAsDouble(reporter, "LSB_stepsize_decrease_factor", stepsize_decrease_factor_);
    options->valueAsString(reporter, "LSB_direction_search_type", direction_search_type_);
    options->valueAsString(reporter, "LSB_directional_derivative_first_order_type",
                           directional_derivative_first_order_type_);
    options->valueAsString(reporter, "LSB_directional_derivative_second_order_type",
                           directional_derivative_second_order_type_);
    options->valueAsString(reporter, "LSB_method_pure_first_order", method_pure_first_order_);
    options->valueAsString(reporter, "LSB_method_pure_second_order", method_pure_second_order_);
    options->valueAsString(reporter, "LSB_method_hybrid", method_hybrid_);
    // options->valueAsInteger(reporter, "linesearch_max_backtrack", max_backtrack_);

}  // end getOptions

// Initialize
void LineSearchBacktracking::initialize(const Options* options, Quantities* quantities, const Reporter* reporter)
{
    std::string partition_strategy;
    options->valueAsString(reporter, "space_partition", partition_strategy);
    if (partition_strategy.compare("FirstOrderPartition") == 0)
    {
        directional_derivative_second_order_type_ = "Not applicable";
        method_pure_second_order_ = "Not applicable";
    };
    quantities->setStepsizeLineSearch(fmax(stepsize_minimum_, stepsize_initial_));
    details_ = "";
    details_ += "|-- search direction type: " + direction_search_type_ + "\n";
    details_ += "|-- directional derivative first order type: " + directional_derivative_first_order_type_ + "\n";
    details_ += "|-- directional derivative second order type: " + directional_derivative_second_order_type_ + "\n";
    if (direction_search_type_.compare("pure") == 0)
    {
        details_ += "|-- Line search method first order direction: " + method_pure_first_order_ + "\n";
        details_ += "|-- Line search method second order direction: " + method_pure_second_order_ + "\n";
    }
    else if (direction_search_type_.compare("hybrid") == 0)
    {
        details_ += "|-- Line search method: " + method_hybrid_ + "\n";
        THROW_EXCEPTION(LS_EVALUATION_FAILURE_EXCEPTION,
                        "No implementation for the LSB_direction_search_type being set to hybrid yet\n");
    }
    else
    {
        details_ += "|-- Line search method: UNKNOWN\n";
        std::string msg = "Invalid option for direction_search_type: " + direction_search_type_ + "\n";
        THROW_EXCEPTION(LS_EVALUATION_FAILURE_EXCEPTION, msg);
    }  // end if
}

// Run line search - get the acceptable trial iterate; move from current iterate to the trial
// iterate is done in the solver
void LineSearchBacktracking::runLineSearch(const Options* options, Quantities* quantities, const Reporter* reporter,
                                           Strategies* strategies)
{
    // Initialize values
    setStatus(LS_UNSET);
    // set number_of_backtrack_ to zero;
    reset();
    // declare some variables for printing purpose
    std::string iteration_type = "desc";
    int         new_zero_groups = 0;
    double      directional_derivative = FARSA_DOUBLE_INFINITY;
    int         projection_attemps = 0;
    std::string search_method = "UNKNOWN";

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

        // determine quantities->direction_search() from the search direction type
        // directional derivative and backtracking linesearch method are also be determined here
        directional_derivative = 0.0;
        auto direction_search = quantities->directionSearch();

        if (direction_search_type_.compare("pure") == 0)
        {
            if ((strategies->directionComputationFirstOrder()->status() == DC_SUCCESS) and
                (strategies->directionComputationSecondOrder()->status() == DC_SKIPPED))
            {
                // compute the search direction
                quantities->directionSearch()->copy(*(quantities->directionFirstOrder()));
                // determine the directional derivative
                if (directional_derivative_first_order_type_.compare("proximal_step") == 0)
                {
                    directional_derivative = direction_search->innerProduct(*direction_search);
                    directional_derivative /= -quantities->stepsizeProximalGradient();
                }
                else
                {
                    std::string msg = "Invalid option for directional_derivative_first_order_type: " +
                                      directional_derivative_first_order_type_ + "\n";
                    THROW_EXCEPTION(LS_EVALUATION_FAILURE_EXCEPTION, msg);
                }
                // determine the search  method
                search_method = method_pure_first_order_;
            }
            else if ((strategies->directionComputationFirstOrder()->status() == DC_SKIPPED) and
                     (strategies->directionComputationSecondOrder()->status() == DC_SUCCESS))
            {
                // determine the directional derivative
                if (directional_derivative_second_order_type_.compare("gradient") == 0)
                {
                    // do not pre-compute the dirder for projected_armijo_groupl1
                    if (method_pure_second_order_.compare("armijo") == 0)
                    {
                        quantities->directionSearch()->copy(*(quantities->directionSecondOrder()));
                        directional_derivative = direction_search->innerProduct(*direction_search);
                    }
                }
                else
                {
                    std::string msg = "Invalid option for directional_derivative_first_order_type: " +
                                      directional_derivative_first_order_type_ + "\n";
                    THROW_EXCEPTION(LS_EVALUATION_FAILURE_EXCEPTION, msg);
                }
                // determine the search  method
                search_method = method_pure_second_order_;
            }
            else
            {
                THROW_EXCEPTION(LS_EVALUATION_FAILURE_EXCEPTION,
                                "Either both directionComputationFirstOrder and directionComputationSecondOrder are "
                                "performed or both are not performed. This is not supposed to happen.\n");
            }
        }
        else if (direction_search_type_.compare("hybrid") == 0)
        {
            // Vector temp(quantities->numberOfVariables());
            // temp.copy(*quantities->directionFirstOrder());
            // temp.addScaledVector(1.0, *quantities->directionSecondOrder());
            // quantities->directionSearch()->copy(temp);
            THROW_EXCEPTION(LS_EVALUATION_FAILURE_EXCEPTION,
                            "No implementation for the LSB_direction_search_type being set to hybrid yet\n");
            search_method = method_hybrid_;
        }
        else
        {
            std::string msg = "Invalid option for direction_search_type: " + direction_search_type_ + "\n";
            THROW_EXCEPTION(LS_EVALUATION_FAILURE_EXCEPTION, msg);
        }  // end if
        // sanity check on the search_method
        if (search_method.compare("UNKNOWN") == 0)
        {
            THROW_EXCEPTION(LS_EVALUATION_FAILURE_EXCEPTION,
                            "The search_method is not specified. This is not supposed to happen.\n");
        }

        // perform the Armijo-backtrack linesearch
        if (search_method.compare("armijo") == 0)
        {
            // Loop
            while (true)
            {
                // Declare new point
                quantities->setTrialIterate(quantities->currentIterate()->makeNewLinearCombination(
                    1.0, quantities->stepsizeLineSearch(), *quantities->directionSearch()));

                // Evaluate trial objective for both smooth and nonsmooth function
                evaluation_success = quantities->trialIterate()->evaluateObjectiveAll(*quantities);

                // Check for successful evaluation
                if (evaluation_success)
                {
                    // Check for sufficient decrease
                    bool sufficient_decrease =
                        (quantities->trialIterate()->objectiveAll() - quantities->currentIterate()->objectiveAll() <=
                         stepsize_sufficient_decrease_threshold_ * quantities->stepsizeLineSearch() *
                             directional_derivative);
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
                        if (quantities->trialIterate()->objectiveAll() < quantities->currentIterate()->objectiveAll())
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
                quantities->setStepsizeLineSearch(stepsize_decrease_factor_ * quantities->stepsizeLineSearch());
                // begin to backtrack
                number_of_backtrack_ += 1;

            }  // end while
        }
        else if (search_method.compare("projected_armijo_groupl1") == 0)
        {
            // determine the intersection of the search direction with the projection ball

            // infity means the corresponding group will never be projected to zero
            std::vector<double> tau_vec(quantities->numberOfGroups(), FARSA_DOUBLE_INFINITY);
            auto   ptr = std::static_pointer_cast<SpacePartitionGroupL1PGBased>(strategies->spacePartition());
            double kappa1 = ptr->kappa1();
            double kappa2 = ptr->kappa2();
            auto   d = quantities->directionSecondOrder();
            auto   x = quantities->currentIterate()->vector();
            double norm_gradient_all_indices_working = quantities->currentIterate()->normGradientAllIndicesWorking();
            auto   per_group_gradient_all_2norm = quantities->currentIterate()->perGroupGradientAll2Norm();
            double temp;
            int    start;
            int    end;
            double norm_dGi_sq;
            double dGiTxGi;
            double xGiTxGi;
            double rho_i;
            double radius_i;
            double a, b, c, delta;
            double tau_min = FARSA_DOUBLE_INFINITY;

            for (auto i : *quantities->groupsWorking())
            {
                start = (*quantities->groups())[i][0];
                end = (*quantities->groups())[i][1];
                norm_dGi_sq = 0.0;
                dGiTxGi = 0.0;
                xGiTxGi = 0.0;
                // calcuate ||d_Gi||^2, d_Gi^T x_Gi, x_Gi^T x_Gi
                for (int j = start; j <= end; j++)
                {
                    norm_dGi_sq += pow(d->values()[j], 2.0);
                    dGiTxGi += d->values()[j] * x->values()[j];
                    xGiTxGi += pow(x->values()[j], 2.0);
                }
                temp = fmin(1.0, kappa2 * norm_gradient_all_indices_working * (end - start + 1) /
                                     quantities->groupsWorking()->size());
                rho_i = fmax(kappa1 * per_group_gradient_all_2norm->values()[i], temp);
                // radius_i = fmin(rho_i, sin(M_PI / 4) * per_group_2norm->values()[i]);
                radius_i = fmin(rho_i, sin(M_PI / 4) * sqrt(xGiTxGi));
                a = norm_dGi_sq;
                b = 2 * dGiTxGi;
                c = xGiTxGi - pow(radius_i, 2.0);
                // compute the first intersection
                delta = pow(b, 2.0) - 4 * a * c;
                if (delta <= 0)
                {
                    tau_vec[i] = FARSA_DOUBLE_INFINITY;
                }
                else
                {
                    temp = fmin((-b + sqrt(delta)) / (2 * a), (-b - sqrt(delta)) / (2 * a));
                    if (temp < 0)
                    {
                        tau_vec[i] = FARSA_DOUBLE_INFINITY;
                    }
                    else
                    {
                        tau_vec[i] = temp;
                        tau_min = fmin(temp, tau_min);
                    }
                }
            }
            // begin backtracking
            projection_attemps = 0;
            new_zero_groups = 0;

            // form the gradient all in the full space and padding zero when necessary
            Vector gradient_all_full_space(quantities->numberOfVariables());
            auto   grad_smooth = quantities->currentIterate()->gradientSmooth();
            auto   grad_nonsmooth = quantities->currentIterate()->gradientNonsmooth();
            for (auto i : *quantities->indiciesWorking())
            {
                gradient_all_full_space.valuesModifiable()[i] = grad_smooth->values()[i] + grad_nonsmooth->values()[i];
            }

            // projection stage: non-increasing
            while (quantities->stepsizeLineSearch() >= tau_min)
            {
                new_zero_groups = 0;
                // Declare new point
                quantities->setTrialIterate(
                    quantities->currentIterate()->makeNewLinearCombination(1.0, 0.0, *quantities->directionSearch()));
                for (int i = 0; i < quantities->numberOfGroups(); i++)
                {
                    start = (*quantities->groups())[i][0];
                    end = (*quantities->groups())[i][1];
                    // current group is projected to 0
                    if (quantities->stepsizeLineSearch() > tau_vec[i])
                    {
                        for (int j = start; j <= end; j++)
                        {
                            quantities->trialIterate()->vector()->valuesModifiable()[j] = 0.0;
                        }
                        new_zero_groups += 1;
                    }
                    else
                    {
                        // take a step along the second order direction
                        for (int j = start; j <= end; j++)
                        {
                            quantities->trialIterate()->vector()->valuesModifiable()[j] =
                                x->values()[j] + quantities->stepsizeLineSearch() * d->values()[j];
                        }
                    }
                }
                projection_attemps += 1;
                // Evaluate trial objective for both smooth and nonsmooth function
                evaluation_success = quantities->trialIterate()->evaluateObjectiveAll(*quantities);

                // Check for successful evaluation
                if (evaluation_success)
                {
                    // Check for sufficient decrease
                    bool non_increasing =
                        (quantities->trialIterate()->objectiveAll() - quantities->currentIterate()->objectiveAll() <=
                         0);

                    // Check Armijo condition
                    if (non_increasing)
                    {
                        iteration_type = "proj";
                        // actual direction taken with unit-stepsize
                        auto search_direction = quantities->trialIterate()->vector()->makeNewLinearCombination(
                            1.0, -1.0, *(quantities->currentIterate()->vector()));
                        quantities->directionSearch()->copy(*(search_direction));
                        quantities->setStepsizeLineSearch(1.0);
                        directional_derivative = search_direction->innerProduct(gradient_all_full_space);
                        THROW_EXCEPTION(LS_SUCCESS_EXCEPTION, "Line search successful.");
                    }  // end if

                }  // end if
                quantities->setStepsizeLineSearch(stepsize_decrease_factor_ * quantities->stepsizeLineSearch());
            }
            // no further progress can be made due to small directional derivative along the second order direciton;
            // terminate!
            quantities->directionSearch()->copy(*d);
            directional_derivative = direction_search->innerProduct(gradient_all_full_space);
            double ratio = fabs(directional_derivative) / (1 + fabs(quantities->currentIterate()->objectiveAll()));
            if (ratio < 1e-15)
            {
                iteration_type = "npgs";
                quantities->setStepsizeLineSearch(1.0);
                THROW_EXCEPTION(LS_NO_FURTHUR_PROGRESS_EXCEPTION,
                                "No further progress can be made. Early termination.");
            }

            // sufficient decrease stage
            while (true)
            {
                // Declare new point
                quantities->setTrialIterate(quantities->currentIterate()->makeNewLinearCombination(
                    1.0, quantities->stepsizeLineSearch(), *quantities->directionSearch()));

                // Evaluate trial objective for both smooth and nonsmooth function
                evaluation_success = quantities->trialIterate()->evaluateObjectiveAll(*quantities);

                // Check for successful evaluation
                if (evaluation_success)
                {
                    // Check for sufficient decrease
                    bool sufficient_decrease =
                        (quantities->trialIterate()->objectiveAll() - quantities->currentIterate()->objectiveAll() <=
                         stepsize_sufficient_decrease_threshold_ * quantities->stepsizeLineSearch() *
                             directional_derivative);
                    // std::cout.precision(17);
                    // std::cout << "\n\nFnew: " << quantities->trialIterate()->objectiveAll()
                    //           << " Fnow: " << quantities->currentIterate()->objectiveAll() << "  const: "
                    //           << stepsize_sufficient_decrease_threshold_ * quantities->stepsizeLineSearch()
                    //           << " dirder: " << directional_derivative << std::endl;
                    // std::cout << "LHS: "
                    //           << quantities->trialIterate()->objectiveAll() -
                    //                  quantities->currentIterate()->objectiveAll()
                    //           << " RHS: "
                    //           << stepsize_sufficient_decrease_threshold_ * quantities->stepsizeLineSearch() *
                    //                  directional_derivative
                    //           << " LHS <= RHS: " << sufficient_decrease << "\n"
                    //           << std::endl;

                    // Check Armijo condition
                    if (sufficient_decrease)
                    {
                        iteration_type = "desc";
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
                        if (quantities->trialIterate()->objectiveAll() < quantities->currentIterate()->objectiveAll())
                        {
                            iteration_type = "small";
                            THROW_EXCEPTION(LS_SUCCESS_EXCEPTION, "Line search successful.");

                        }  // end if

                    }  // end if

                    // Set null step
                    quantities->setStepsizeLineSearch(0.0);

                    // Set new point
                    quantities->setTrialIterateToCurrentIterate();

                    // Terminate
                    iteration_type = "small";
                    THROW_EXCEPTION(LS_SUCCESS_EXCEPTION, "Line search successful.");

                }  // end if

                // Update stepsize
                quantities->setStepsizeLineSearch(stepsize_decrease_factor_ * quantities->stepsizeLineSearch());
                // begin to backtrack
                number_of_backtrack_ += 1;

            }  // end while
        }
        else if (search_method.compare("dogleg") == 0)
        {
            THROW_EXCEPTION(LS_EVALUATION_FAILURE_EXCEPTION,
                            "No implementation for the search method being set to "
                            "dogleg yet\n");
        }
        else
        {
            std::string msg = "Invalid option for search method: " + search_method + "\n";
            THROW_EXCEPTION(LS_EVALUATION_FAILURE_EXCEPTION, msg);
        }

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
    }
    catch (LS_NO_FURTHUR_PROGRESS_EXCEPTION& exec)
    {
        setStatus(LS_NO_FURTHUR_PROGRESS);
    }
    // end catch

    if (verbose_)
    {
        if (strategies->spacePartition()->name().compare("GroupL1PGBasedPartition") == 0)
        {
            if (search_method.compare("projected_armijo_groupl1") == 0)
            {
                reporter->printf(R_SOLVER, R_PER_ITERATION, " %+.2e  %s %6d %+.2e %2d/%2d",
                                 quantities->directionSecondOrder()->norm2(), iteration_type.c_str(), new_zero_groups,
                                 directional_derivative, projection_attemps, number_of_backtrack_);
            }
            else
            {
                reporter->printf(R_SOLVER, R_PER_ITERATION, " %+.2e  %s %s %+.2e %s/%2d",
                                 quantities->directionSearch()->norm2(), iteration_type.c_str(), "------",
                                 directional_derivative, "--", number_of_backtrack_);
            }
        }
        else
        {
            if (search_method.compare("projected_armijo_groupl1"))
            {
                reporter->printf(R_SOLVER, R_PER_ITERATION, " %+.2e", quantities->directionSecondOrder()->norm2());
            }
            else
            {
                reporter->printf(R_SOLVER, R_PER_ITERATION, " %+.2e", quantities->directionSearch()->norm2());
            }
        }
    }
    else
    {
        // Print iteration information
        reporter->printf(R_SOLVER, R_PER_ITERATION, " %+.2e", quantities->directionSearch()->norm2());
    }
    reporter->printf(R_SOLVER, R_PER_ITERATION, "  %+.2e", quantities->stepsizeLineSearch());

}  // end runLineSearch

}  // namespace FaRSA
