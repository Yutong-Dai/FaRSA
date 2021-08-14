// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __TESTSOLVERSECONDORDER_HPP__
#define __TESTSOLVERSECONDORDER_HPP__

#include <iostream>

#include "FaRSAGroupL1.hpp"
#include "FaRSALinearRegressionLoss.hpp"
#include "FaRSASolver.hpp"
using namespace FaRSA;

// Implementation of test
int testSolverSecondOrderImplementation(int option)
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
        new LinearRegressionLoss((char*)"./data/python_impl/FaRSAGroup/test/lsMatrix.txt", M_COORDINATE_LIST,
                                 (char*)"./data/python_impl/FaRSAGroup/test/lsLabel.txt", "bodyfat"));

    // create non-smooth funciton
    // set-up groups
    std::shared_ptr<std::vector<std::vector<int>>> groups(new std::vector<std::vector<int>>);
    std::vector<int>                               temp;
    temp.push_back(0);
    temp.push_back(5);
    groups->push_back(temp);
    temp.clear();
    temp.push_back(6);
    temp.push_back(7);
    groups->push_back(temp);
    temp.clear();
    temp.push_back(8);
    temp.push_back(12);
    groups->push_back(temp);
    temp.clear();
    temp.push_back(13);
    temp.push_back(13);
    groups->push_back(temp);
    std::shared_ptr<std::vector<double>> weights(new std::vector<double>);
    weights->push_back(sqrt(6.0));
    weights->push_back(sqrt(2.0));
    weights->push_back(sqrt(5.0));
    weights->push_back(sqrt(1.0));
    double                   penaty = 0.1;
    std::shared_ptr<GroupL1> r(new GroupL1(groups, weights, penaty));

    // create a starting point
    std::shared_ptr<Vector> x(new Vector(6));
    // tests begins here
    // Vector ans;
    // Declare solver object
    FaRSASolver farsa;

    // Modify options from file
    farsa.options()->modifyOptionsFromFile(&reporter, "farsa.opt");

    // Optimize
    farsa.optimize(f, r, *x);

    // // Check option
    // if (option == 1)
    // {
    //     // Print final message
    //     if (result == 0)
    //     {
    //         reporter.printf(R_SOLVER, R_BASIC, "TEST WAS SUCCESSFUL.\n");
    //     }
    //     else
    //     {
    //         reporter.printf(R_SOLVER, R_BASIC, "TEST FAILED.\n");
    //     }

    // }  // end if
    // Return
    return result;

}  // end testSolverSecondOrderImplementation

#endif /* __TESTSOLVERSECONDORDER_HPP__ */
