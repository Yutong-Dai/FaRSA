// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __TESTLINEARREGRESSIONLOSS_HPP__
#define __TESTLINEARREGRESSIONLOSS_HPP__

#include <algorithm>
#include <cmath>
#include <iostream>

#include "FaRSALinearRegressionLoss.hpp"
#include "FaRSAReporter.hpp"
using namespace FaRSA;

// Implementation of test
int testLinearRegressionLossImplementation(int option)
{
    // Initialize output
    int result = 0;

    // Declare reporter
    Reporter reporter;

    // Check option
    if (option == 1)
    {
        // Declare stream report
        std::shared_ptr<StreamReport> s(
            new StreamReport("s", R_SOLVER, R_BASIC));

        // Set stream report to standard output
        s->setStream(&std::cout);

        // Add stream report to reporter
        reporter.addReport(s);

    }  // end if

    // // Declare matrix
    // Matrix A;
    // Vector b;
    // char*  file_name;
    // //  Read from file name
    // file_name = (char*)"./data/lsMatrix.txt";
    // A.setFromFile(file_name, M_COORDINATE_LIST);
    // A.print(&reporter, "Matrix A read from file:");
    // // Print matrix
    // file_name = (char*)"./data/lsLabel.txt";
    // b.setFromFile(file_name);
    // b.print(&reporter, "Vector b read from file:");

    // create  LinearRegressionLoss
    LinearRegressionLoss ls((char*)"./data/lsMatrix.txt", M_COORDINATE_LIST,
                            (char*)"./data/lsLabel.txt", "stdNormal106");
    double               f;
    Vector               x(6, 4.17);
    Vector               g(6);
    std::vector<int>     cols({0, 5, 4, 2});
    Vector               Hv(4);
    Vector               v(4);
    std::vector<float>   gans({5.98971467, 0.43315364, 3.4363911, 0.56336833,
                             0.80778718, 0.20332681});
    std::vector<float>   Hvans({0.5799699, 1.22246231, 0.43032134, 0.13205128});

    for (int i = 0; i < v.length(); i++)
    {
        v.valuesModifiable()[i] = i * 1.1;
    }
    // test function value evaluation
    bool fevalSuccess = ls.evaluateObjective(x, f);
    if (!fevalSuccess or fabs(f - 24.163181135723303) > 1e-7)
    {
        result = 1;
        reporter.printf(R_SOLVER, R_BASIC,
                        "Expected fval:%8.5f | Actual fval:%8.5f\n",
                        24.163181135723303, f);
    }

    // test gradient evaluation
    std::vector<int> placeholder;
    bool gevalSuccess = ls.evaluateGradient(x, placeholder, g);
    if (!gevalSuccess)
    {
        result = 1;
    }
    for (int i = 0; i < g.length(); i++)
    {
        if (fabs(g.values()[i] - gans[i]) > 1e-7)
        {
            reporter.printf(R_SOLVER, R_BASIC,
                            "Expected g[%d]:%8.5f | Actual g[%d]:%8.5f | "
                            "Difference: %20.16f\n",
                            gans[i], i, g.values()[i], i,
                            fabs(g.values()[i] - gans[i]));
            result = 1;
        }
    }

    ls.evaluateGradient(x, cols, g);

    // test Hessian-Vector Product for submatrix
    std::sort(cols.begin(), cols.end());
    bool hvevalSuccess = ls.evaluateHessianVectorProduct(x, cols, v, Hv);

    if (!hvevalSuccess)
    {
        result = 1;
    }
    for (int i = 0; i < Hv.length(); i++)
    {
        if (fabs(Hv.values()[i] - Hvans[i]) > 1e-7)
        {
            reporter.printf(R_SOLVER, R_BASIC,
                            "Expected Hv[%d]:%8.5f | Actual Hv[%d]:%8.5f | "
                            "Difference: %20.16f\n",
                            Hvans[i], i, Hv.values()[i], i,
                            fabs(Hv.values()[i] - Hvans[i]));
            result = 1;
        }
    }
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

}  // end testMatrixImplementation

#endif /* __TESTLINEARREGRESSIONLOSS_HPP__ */
