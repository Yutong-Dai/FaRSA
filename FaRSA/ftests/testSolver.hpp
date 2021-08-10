// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __TESTSOLVER_HPP__
#define __TESTSOLVER_HPP__

#include <iostream>

#include "FaRSAGroupL1.hpp"
#include "FaRSALinearRegressionLoss.hpp"
#include "FaRSASolver.hpp"
using namespace FaRSA;

// Implementation of test
int testSolverImplementation(int option)
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

    // create a starting point
    std::shared_ptr<Vector> x(new Vector(6));
    x->valuesModifiable()[2] = 1.1;
    x->valuesModifiable()[3] = 0.0;
    x->valuesModifiable()[4] = 2.2;
    x->valuesModifiable()[5] = 3.3;
    // tests begins here
    // Vector ans;
    // Declare solver object
    FaRSASolver farsa;

    // Modify options from file
    // farsa.options()->modifyOptionsFromFile("nonopt.opt");

    // Optimize
    farsa.optimize(f, r, *x);

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

}  // end testSolverImplementation

#endif /* __TESTSOLVER_HPP__ */
