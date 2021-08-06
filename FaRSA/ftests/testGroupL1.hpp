// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __TESTGROUPL1_HPP__
#define __TESTGROUPL1_HPP__

#include <FaRSAQuantities.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>

#include "FaRSAGroupL1.hpp"
#include "FaRSAReporter.hpp"
using namespace FaRSA;

// Implementation of test
int testGroupL1Implementation(int option)
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

    // set-up groups
    std::shared_ptr<std::vector<std::vector<int>>> groups(
        new std::vector<std::vector<int>>);
    std::vector<int> temp;
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
    // std::vector<std::vector<int>> groups{{0, 1}, {2, 4}, {5, 5}};

    // std::vector<double> weights;
    // weights.push_back(sqrt(2.0));
    // weights.push_back(sqrt(3.0));
    // weights.push_back(sqrt(1.0));
    std::shared_ptr<std::vector<double>> weights(new std::vector<double>);
    weights->push_back(sqrt(2.0));
    weights->push_back(sqrt(3.0));
    weights->push_back(sqrt(1.0));
    double penaty = 1.23;

    GroupL1 r(groups, weights, penaty);
    Vector  ans;

    // test function value evaluation
    double f;
    Vector x(6);

    x.valuesModifiable()[2] = 1.1;
    x.valuesModifiable()[3] = 0.0;
    x.valuesModifiable()[4] = 2.2;
    x.valuesModifiable()[5] = 3.3;
    bool fevalSuccess = r.evaluateObjective(x, f);
    // reporter.printf(R_SOLVER, R_BASIC, "func %e\n", f);
    ans.setFromFile((char*)"./data/rfun.txt");
    if (!fevalSuccess or fabs(f - ans.values()[0]) > 1e-7)
    {
        result = 1;
        reporter.printf(R_SOLVER, R_BASIC,
                        "Expected fval:%8.5f | Actual fval:%8.5f\n",
                        ans.values()[0], f);
    }

    // test gradient evaluation
    std::vector<int> cols({2, 3, 4, 5});
    Vector           g(cols.size());
    bool             gevalSuccess = r.evaluateGradient(x, cols, g);
    ans.setFromFile((char*)"./data/rgrad.txt");
    // g.print(&reporter, "gradient:");
    if (!gevalSuccess)
    {
        result = 1;
    }
    for (int i = 0; i < g.length(); i++)
    {
        if (fabs(g.values()[i] - ans.values()[i]) > 1e-7)
        {
            reporter.printf(R_SOLVER, R_BASIC,
                            "Expected g[%d]:%8.5f | Actual g[%d]:%8.5f | "
                            "Difference: %20.16f\n",
                            ans.values()[i], i, g.values()[i], i,
                            fabs(g.values()[i] - ans.values()[i]));
            result = 1;
        }
    }

    // test proximal gradient evaluation
    Quantities quantities;
    quantities.setStepsize(0.2);
    Vector gradfx(6);
    Vector proxgrad(6);
    for (int i = 0; i < gradfx.length(); i++)
    {
        gradfx.valuesModifiable()[i] = x.values()[i] + 0.1;
    }
    // gradfx.print(&reporter, "gradfx:");
    bool proxevalSuccess;
    proxevalSuccess =
        r.computeProximalGradientUpdate(x, gradfx, quantities, proxgrad);
    ans.setFromFile((char*)"./data/rprox1.txt");
    // proxgrad.print(&reporter, "proximal-grad:");
    if (!proxevalSuccess)
    {
        result = 1;
    }
    for (int i = 0; i < proxgrad.length(); i++)
    {
        if (fabs(proxgrad.values()[i] - ans.values()[i]) > 1e-7)
        {
            reporter.printf(R_SOLVER, R_BASIC,
                            "Expected prox[%d]:%8.5f | Actual prox[%d]:%8.5f | "
                            "Difference: %20.16f\n",
                            ans.values()[i], i, proxgrad.values()[i], i,
                            fabs(proxgrad.values()[i] - ans.values()[i]));
            result = 1;
        }
    }
    for (int i = 0; i < gradfx.length(); i++)
    {
        gradfx.valuesModifiable()[i] = x.values()[i] + 2.0;
    }

    proxevalSuccess =
        r.computeProximalGradientUpdate(x, gradfx, quantities, proxgrad);
    ans.setFromFile((char*)"./data/rprox2.txt");
    if (!proxevalSuccess)
    {
        result = 1;
    }
    for (int i = 0; i < proxgrad.length(); i++)
    {
        if (fabs(proxgrad.values()[i] - ans.values()[i]) > 1e-7)
        {
            reporter.printf(R_SOLVER, R_BASIC,
                            "Expected prox[%d]:%8.5f | Actual prox[%d]:%8.5f | "
                            "Difference: %20.16f\n",
                            ans.values()[i], i, proxgrad.values()[i], i,
                            fabs(proxgrad.values()[i] - ans.values()[i]));
            result = 1;
        }
    }

    // test Hessian-Vector Product for submatrix
    Vector Hv(4);
    Vector v(4);
    for (int i = 0; i < v.length(); i++)
    {
        v.valuesModifiable()[i] = i * 1.1;
    }
    std::sort(cols.begin(), cols.end());
    bool hvevalSuccess = r.evaluateHessianVectorProduct(x, cols, v, Hv);
    ans.setFromFile((char*)"./data/rHv.txt");
    if (!hvevalSuccess)
    {
        result = 1;
    }
    for (int i = 0; i < Hv.length(); i++)
    {
        if (fabs(Hv.values()[i] - ans.values()[i]) > 1e-7)
        {
            reporter.printf(R_SOLVER, R_BASIC,
                            "Expected Hv[%d]:%8.5f | Actual Hv[%d]:%8.5f | "
                            "Difference: %20.16f\n",
                            ans.values()[i], i, proxgrad.values()[i], i,
                            fabs(Hv.values()[i] - ans.values()[i]));
            result = 1;
        }
    }
    // Hv.print(&reporter, "Hv:");
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

}  // end testGroupL1Implementation

#endif /* __TESTGROUPL1_HPP__ */
