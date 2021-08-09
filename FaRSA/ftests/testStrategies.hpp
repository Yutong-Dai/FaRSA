// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __TESTSTRATIGIES_HPP__
#define __TESTSTRATIGIES_HPP__

#include <FaRSAQuantities.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>

#include "FaRSADirectionComputationProximalGradient.hpp"
#include "FaRSAGroupL1.hpp"
#include "FaRSALineSearchBacktracking.hpp"
#include "FaRSALinearRegressionLoss.hpp"
#include "FaRSAOptions.hpp"
#include "FaRSAPoint.hpp"
#include "FaRSAReporter.hpp"
#include "FaRSASpacePartitionFirstOrder.hpp"
using namespace FaRSA;

// Implementation of test
int testStrategiesImplementation(int option)
{
    // Initialize output
    int result = 0;

    // Declare reporter
    Reporter reporter;

    // Check option
    if (option == 1)
    {
        // Declare stream report
        std::shared_ptr<StreamReport> s(new StreamReport("s", R_SOLVER, R_BASIC));

        // Set stream report to standard output
        s->setStream(&std::cout);

        // Add stream report to reporter
        reporter.addReport(s);

    }  // end if

    // create smooth function
    std::shared_ptr<LinearRegressionLoss> f(
        new LinearRegressionLoss((char*)"./data/lsMatrix.txt", M_COORDINATE_LIST,
                                 (char*)"./data/lsLabel.txt", "stdNormal106"));

    // create non-smooth funciton
    // set-up groups
    std::shared_ptr<std::vector<std::vector<int>>> groups(new std::vector<std::vector<int>>);
    std::vector<int>                               temp;
    temp.push_back(0);
    temp.push_back(1);
    groups->push_back(temp);
    temp.clear();
    temp.push_back(2);
    temp.push_back(4);
    groups->push_back(temp);
    temp.clear();
    temp.push_back(5);
    temp.push_back(5);
    groups->push_back(temp);
    std::shared_ptr<std::vector<double>> weights(new std::vector<double>);
    weights->push_back(sqrt(2.0));
    weights->push_back(sqrt(3.0));
    weights->push_back(sqrt(1.0));
    double                   penaty = 1.23;
    std::shared_ptr<GroupL1> r(new GroupL1(groups, weights, penaty));

    // Declare point
    std::shared_ptr<Vector> x(new Vector(6));
    x->valuesModifiable()[2] = 1.1;
    x->valuesModifiable()[3] = 0.0;
    x->valuesModifiable()[4] = 2.2;
    x->valuesModifiable()[5] = 3.3;
    Point      p(f, r, x, 1.0);
    Quantities quantities;
    Options    options;
    Strategies strategies;
    quantities.addOptions(&options, &reporter);
    quantities.getOptions(&options, &reporter);

    strategies.addOptions(&options, &reporter);
    strategies.getOptions(&options, &reporter);
    // tests begins here
    Vector ans;
    bool   operation_success;

    // print all options and strategies
    // options.print(&reporter);
    strategies.printHeader(&reporter);

    // check initialization for quantities
    operation_success = quantities.initialize(f, r, *x);
    quantities.setScalingThreshold(200.0);
    if (!operation_success)
    {
        result = 1;
        reporter.printf(R_SOLVER, R_BASIC, "Quantities initialization failed.\n");
    }

    // check some member values

    // check initial stepsize
    if (fabs(quantities.stepsizeLineSearch() - 0.0) > 1e-16)
    {
        result = 1;
        reporter.printf(
            R_SOLVER, R_BASIC,
            "The stepsize for linesearch should have been initialized to 0.0, but wasn't.\n");
    }

    if (fabs(quantities.stepsizeProximalGradient() - 1.0) > 1e-16)
    {
        result = 1;
        reporter.printf(R_SOLVER, R_BASIC,
                        "The stepsize for proximal gradient update should have been initialized to "
                        "1.0, but wasn't.\n");
    }
    reporter.printf(R_SOLVER, R_BASIC, "The scale used for the problem is %3.8e.\n",
                    quantities.scaleApplied());

    // compute space partition
    strategies.initialize(&options, &quantities, &reporter);
    strategies.spacePartition()->partitionSpace(&options, &quantities, &reporter, &strategies);
    for (auto element : *quantities.groupsFirstOrder())
    {
        reporter.printf(R_SOLVER, R_BASIC, "groups in first order set: %d\n", element);
    }
    for (auto element : *quantities.indiciesWorking())
    {
        reporter.printf(R_SOLVER, R_BASIC, "working indices: %d\n", element);
    }

    // compute search direction
    quantities.setStepsizeProximalGradient(0.2);
    strategies.directionComputationFirstOrder()->computeDirection(&options, &quantities, &reporter,
                                                                  &strategies);
    reporter.printf(
        R_SOLVER, R_BASIC, "perform directionComputationFirstOrder: %s\n",
        strategies.directionComputationFirstOrder()->performCompuation() ? "true" : "false");
    strategies.directionComputationSecondOrder()->computeDirection(&options, &quantities, &reporter,
                                                                   &strategies);
    reporter.printf(
        R_SOLVER, R_BASIC, "perform directionComputationSecondOrder: %s\n",
        strategies.directionComputationSecondOrder()->performCompuation() ? "true" : "false");
    ans.setFromFile((char*)"./data/trueProx.txt");
    ans.addScaledVector(-1.0, *(quantities.currentIterate()->vector()));
    auto prox = quantities.direction();
    for (int i = 0; i < prox->length(); i++)
    {
        if (fabs(prox->values()[i] - ans.values()[i]) > 1e-7)
        {
            reporter.printf(R_SOLVER, R_BASIC,
                            "Expected prox[%d]:%8.5f | Actual prox[%d]:%8.5f | "
                            "Difference: %20.16f\n",
                            ans.values()[i], i, prox->values()[i], i,
                            fabs(prox->values()[i] - ans.values()[i]));
            result = 1;
        }
    }

    // perform linesearch
    strategies.lineSearch()->runLineSearch(&options, &quantities, &reporter, &strategies);
    auto ptr = std::static_pointer_cast<LineSearchBacktracking>(strategies.lineSearch());
    reporter.printf(R_SOLVER, R_BASIC, "number of backtracks: %d\n", ptr->numberOfBacktrack());
    ans.setFromFile((char*)"./data/backtrack_trial_iterates.txt");
    for (int i = 0; i < ans.length(); i++)
    {
        if (fabs(quantities.trialIterate()->vector()->values()[i] - ans.values()[i]) > 1e-7)
        {
            reporter.printf(
                R_SOLVER, R_BASIC,
                "Expected trialIterate[%d]:%8.5f | Actual trialIterate[%d]:%8.5f | "
                "Difference: %20.16f\n",
                ans.values()[i], i, quantities.trialIterate()->vector()->values()[i], i,
                fabs(quantities.trialIterate()->vector()->values()[i] - ans.values()[i]));
            result = 1;
        }
    }

    // perform parameter update
    auto update_rule_pg = strategies.parameterUpdatePGStepsize();
    reporter.printf(R_SOLVER, R_BASIC, "parameterUpdatePGStepsize: %s\n",
                    update_rule_pg->method().c_str());
    reporter.printf(R_SOLVER, R_BASIC, "stepsize pg before update: %3.3e",
                    quantities.stepsizeProximalGradient());
    update_rule_pg->update(&options, &quantities, &reporter, &strategies);
    reporter.printf(R_SOLVER, R_BASIC, "stepsize pg after update: %3.3e\n",
                    quantities.stepsizeProximalGradient());
    // reset stepsize proximal gradient
    quantities.setStepsizeProximalGradient(0.2);
    update_rule_pg->setMethod("Heuristic");
    reporter.printf(R_SOLVER, R_BASIC, "parameterUpdatePGStepsize: %s\n",
                    update_rule_pg->method().c_str());
    reporter.printf(R_SOLVER, R_BASIC, "stepsize pg before update: %3.3e",
                    quantities.stepsizeProximalGradient());
    update_rule_pg->update(&options, &quantities, &reporter, &strategies);
    reporter.printf(R_SOLVER, R_BASIC, "stepsize pg after update: %3.3e\n",
                    quantities.stepsizeProximalGradient());
    // reset stepsize proximal gradient
    quantities.setStepsizeProximalGradient(0.2);
    update_rule_pg->setMethod("LipschitzEstimate");
    reporter.printf(R_SOLVER, R_BASIC, "parameterUpdatePGStepsize: %s\n",
                    update_rule_pg->method().c_str());
    reporter.printf(R_SOLVER, R_BASIC, "stepsize pg before update: %3.3e",
                    quantities.stepsizeProximalGradient());
    update_rule_pg->update(&options, &quantities, &reporter, &strategies);
    reporter.printf(R_SOLVER, R_BASIC, "stepsize pg after update: %3.3e\n",
                    quantities.stepsizeProximalGradient());
    // Check option
    if (option == 1)
    {
        // Print final message
        if (result == 0)
        {
            reporter.printf(R_SOLVER, R_BASIC, "TEST WAS SUCCESSFUL.\n");
        }
        else
        {
            reporter.printf(R_SOLVER, R_BASIC, "TEST FAILED.\n");
        }

    }  // end if
    // Return
    return result;

}  // end testStrategiesImplementation

#endif /* __TESTSTRATIGIES_HPP__ */
