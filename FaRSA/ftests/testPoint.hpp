// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __TESTPOINT_HPP__
#define __TESTPOINT_HPP__

#include <FaRSAQuantities.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>

#include "FaRSAGroupL1.hpp"
#include "FaRSALinearRegressionLoss.hpp"
#include "FaRSAOptions.hpp"
#include "FaRSAPoint.hpp"
#include "FaRSAReporter.hpp"
using namespace FaRSA;

// Implementation of test
int testPointImplementation(int option)
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

    // create smooth function
    std::shared_ptr<LinearRegressionLoss> f(new LinearRegressionLoss(
        (char*)"./data/lsMatrix.txt", M_COORDINATE_LIST,
        (char*)"./data/lsLabel.txt", "stdNormal106"));

    // create non-smooth funciton
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
    quantities.addOptions(&options, &reporter);
    quantities.setStepsize(0.2);
    quantities.getOptions(&options, &reporter);
    quantities.setScalingThreshold(200.0);
    std::shared_ptr<std::vector<int>> indicies_working(new std::vector<int>);
    indicies_working->push_back(2);
    indicies_working->push_back(3);
    indicies_working->push_back(4);
    indicies_working->push_back(5);
    quantities.setIndiciesWorking(indicies_working);
    // tests begins here
    Vector ans;

    reporter.printf(R_SOLVER, R_BASIC, "Test scale=1.0 computation...\n");
    // determin the scale based on the point p
    bool fSmoothevalSuccess;
    bool gSmoothevalSuccess;
    bool fNonsmoothevalSuccess;
    bool gNonsmoothevalSuccess;
    fSmoothevalSuccess = p.evaluateObjectiveSmooth(quantities);
    gSmoothevalSuccess = p.evaluateGradientSmooth(quantities);
    p.determineScale(quantities);
    fNonsmoothevalSuccess = p.evaluateObjectiveNonsmooth(quantities);
    gNonsmoothevalSuccess = p.evaluateGradientNonsmooth(quantities);

    reporter.printf(R_SOLVER, R_BASIC,
                    "Scale used based on the gradient of the smooth function "
                    "at the initial point is: %9.8f\n",
                    p.scale());

    // check values
    ans.setFromFile((char*)"./data/ffun.txt");
    if (!fSmoothevalSuccess or
        fabs(p.objectiveSmoothUnscaled() - ans.values()[0]) > 1e-7)
    {
        result = 1;
        reporter.printf(R_SOLVER, R_BASIC,
                        "Test Smooth Objective: Expected fval:%8.5f | Actual "
                        "fval:%8.5f | Difference: %8.5f\n",
                        ans.values()[0], p.objectiveSmoothUnscaled(),
                        fabs(p.objectiveSmoothUnscaled() - ans.values()[0]));
    }

    ans.setFromFile((char*)"./data/rfun.txt");
    if (!fNonsmoothevalSuccess or
        fabs(p.objectiveNonsmoothUnscaled() - ans.values()[0]) > 1e-7)
    {
        result = 1;
        reporter.printf(
            R_SOLVER, R_BASIC,
            "Test Nonsmooth Objective: Expected fval:%8.5f | Actual "
            "fval:%8.5f | Difference: %8.5f\n",
            ans.values()[0], p.objectiveNonsmoothUnscaled(),
            fabs(p.objectiveNonsmoothUnscaled() - ans.values()[0]));
    }

    // gradient smooth
    ans.setFromFile((char*)"./data/fgrad.txt");
    auto g = p.gradientSmooth();
    if ((!gSmoothevalSuccess) || (g->length() != x->length()))
    {
        result = 1;
    }

    for (int i = 0; i < g->length(); i++)
    {
        if (fabs((*g).values()[i] - ans.values()[i]) > 1e-7)
        {
            reporter.printf(R_SOLVER, R_BASIC,
                            "Test Smooth gradient: Expected g[%d]:%8.5f | "
                            "Actual g[%d]:%8.5f | "
                            "Difference: %20.16f\n",
                            ans.values()[i], i, (*g).values()[i], i,
                            fabs((*g).values()[i] - ans.values()[i]));
            result = 1;
        }
    }
    // g->print(&reporter, "gf:");
    // gradient nonsmooth
    ans.setFromFile((char*)"./data/rgrad.txt");
    auto gr = p.gradientNonsmooth();
    if ((!gNonsmoothevalSuccess) ||
        (gr->length() != quantities.indiciesWorking()->size()))
    {
        result = 1;
    }

    // gr->print(&reporter, "gr:");
    for (int i = 0; i < gr->length(); i++)
    {
        if (fabs((*gr).values()[i] - ans.values()[i]) > 1e-7)
        {
            reporter.printf(R_SOLVER, R_BASIC,
                            "Test Nonsmooth gradient: Expected g[%d]:%8.5f | "
                            "Actual g[%d]:%8.5f | "
                            "Difference: %20.16f\n",
                            ans.values()[i], i, (*gr).values()[i], i,
                            fabs((*gr).values()[i] - ans.values()[i]));
            result = 1;
        }
    }

    // test Hessian-Vector Product
    // reporter.printf(R_SOLVER, R_BASIC, "hess vector limit: %d | counter:
    // %d\n",
    //                 quantities.hessianVectorProductEvaluationLimit(),
    //                 quantities.hessianVectorCounter());

    std::shared_ptr<Vector> v(new Vector(4, 0.0));
    for (int i = 0; i < v->length(); i++)
    {
        v->valuesModifiable()[i] = i * 1.1;
    }
    bool hvSmoothEvalSuccess;
    hvSmoothEvalSuccess = p.evaluateHessianVectorProductSmooth(v, quantities);
    ans.setFromFile((char*)"./data/fHv.txt");
    auto fHv = p.hessianVectorProductSmooth();
    if ((!hvSmoothEvalSuccess) || (fHv->length() != v->length()))
    {
        result = 1;
    }
    for (int i = 0; i < fHv->length(); i++)
    {
        if (fabs(fHv->values()[i] - ans.values()[i]) > 1e-7)
        {
            reporter.printf(R_SOLVER, R_BASIC,
                            "Expected Hv[%d]:%8.5f | Actual Hv[%d]:%8.5f | "
                            "Difference: %20.16f\n",
                            ans.values()[i], i, fHv->values()[i], i,
                            fabs(fHv->values()[i] - ans.values()[i]));
            result = 1;
        }
    }

    bool hvNonsmoothEvalSuccess;
    hvNonsmoothEvalSuccess =
        p.evaluateHessianVectorProductNonsmooth(v, quantities);
    ans.setFromFile((char*)"./data/rHv.txt");
    auto rHv = p.hessianVectorProductNonsmooth();
    if ((!hvNonsmoothEvalSuccess) || (rHv->length() != v->length()))
    {
        result = 1;
    }
    for (int i = 0; i < rHv->length(); i++)
    {
        if (fabs(rHv->values()[i] - ans.values()[i]) > 1e-7)
        {
            reporter.printf(R_SOLVER, R_BASIC,
                            "Expected Hv[%d]:%8.5f | Actual Hv[%d]:%8.5f | "
                            "Difference: %20.16f\n",
                            ans.values()[i], i, rHv->values()[i], i,
                            fabs(rHv->values()[i] - ans.values()[i]));
            result = 1;
        }
    }

    // test proximal gradient
    bool proxComputeSuccess;
    proxComputeSuccess = p.computeProximalGradientUpdate(quantities);
    ans.setFromFile((char*)"./data/trueProx.txt");
    auto prox = p.proximalGraidentUpdate();
    if ((!proxComputeSuccess) || (prox->length() != x->length()))
    {
        result = 1;
    }
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
    reporter.printf(R_SOLVER, R_BASIC, "Test scaled computation...\n");
    Point pnew(f, r, x, 1.0);
    quantities.setScalingThreshold(0.1);
    fSmoothevalSuccess = pnew.evaluateObjectiveSmooth(quantities);
    gSmoothevalSuccess = pnew.evaluateGradientSmooth(quantities);
    pnew.determineScale(quantities);
    fNonsmoothevalSuccess = pnew.evaluateObjectiveNonsmooth(quantities);
    gNonsmoothevalSuccess = pnew.evaluateGradientNonsmooth(quantities);
    reporter.printf(R_SOLVER, R_BASIC,
                    "Scale used based on the gradient of the smooth function "
                    "at the initial point is: %9.8f\n",
                    pnew.scale());

    // check values
    ans.setFromFile((char*)"./data/ffun_scaled.txt");
    if (!fSmoothevalSuccess or
        fabs(pnew.objectiveSmooth() - ans.values()[0]) > 1e-7)
    {
        result = 1;
        reporter.printf(R_SOLVER, R_BASIC,
                        "Test Smooth Objective: Expected fval:%8.5f | Actual "
                        "fval:%8.5f | Difference: %8.5f\n",
                        ans.values()[0], pnew.objectiveSmooth(),
                        fabs(pnew.objectiveSmoothUnscaled() - ans.values()[0]));
    }

    ans.setFromFile((char*)"./data/rfun_scaled.txt");
    if (!fNonsmoothevalSuccess or
        fabs(pnew.objectiveNonsmooth() - ans.values()[0]) > 1e-7)
    {
        result = 1;
        reporter.printf(
            R_SOLVER, R_BASIC,
            "Test Nonsmooth Objective: Expected fval:%8.5f | Actual "
            "fval:%8.5f | Difference: %8.5f\n",
            ans.values()[0], pnew.objectiveNonsmooth(),
            fabs(pnew.objectiveNonsmooth() - ans.values()[0]));
    }

    // gradient smooth
    ans.setFromFile((char*)"./data/fgrad_scaled.txt");
    auto gnew = pnew.gradientSmooth();
    if ((!gSmoothevalSuccess) || (gnew->length() != x->length()))
    {
        result = 1;
    }

    for (int i = 0; i < gnew->length(); i++)
    {
        if (fabs((*gnew).values()[i] - ans.values()[i]) > 1e-7)
        {
            reporter.printf(R_SOLVER, R_BASIC,
                            "Test Smooth gradient: Expected g[%d]:%8.5f | "
                            "Actual g[%d]:%8.5f | "
                            "Difference: %20.16f\n",
                            ans.values()[i], i, (*gnew).values()[i], i,
                            fabs((*gnew).values()[i] - ans.values()[i]));
            result = 1;
        }
    }
    // gradient nonsmooth
    ans.setFromFile((char*)"./data/rgrad_scaled.txt");
    auto grnew = pnew.gradientNonsmooth();
    if ((!gNonsmoothevalSuccess) ||
        (grnew->length() != quantities.indiciesWorking()->size()))
    {
        result = 1;
    }

    for (int i = 0; i < grnew->length(); i++)
    {
        if (fabs((*grnew).values()[i] - ans.values()[i]) > 1e-7)
        {
            reporter.printf(R_SOLVER, R_BASIC,
                            "Test Nonsmooth gradient: Expected g[%d]:%8.5f | "
                            "Actual g[%d]:%8.5f | "
                            "Difference: %20.16f\n",
                            ans.values()[i], i, (*grnew).values()[i], i,
                            fabs((*grnew).values()[i] - ans.values()[i]));
            result = 1;
        }
    }

    // test Hessian-Vector Product
    hvSmoothEvalSuccess =
        pnew.evaluateHessianVectorProductSmooth(v, quantities);
    ans.setFromFile((char*)"./data/fHv_scaled.txt");
    auto fHvnew = pnew.hessianVectorProductSmooth();
    if ((!hvSmoothEvalSuccess) || (fHvnew->length() != v->length()))
    {
        result = 1;
    }
    for (int i = 0; i < fHvnew->length(); i++)
    {
        if (fabs(fHvnew->values()[i] - ans.values()[i]) > 1e-7)
        {
            reporter.printf(R_SOLVER, R_BASIC,
                            "Expected Hv[%d]:%8.5f | Actual Hv[%d]:%8.5f | "
                            "Difference: %20.16f\n",
                            ans.values()[i], i, fHvnew->values()[i], i,
                            fabs(fHvnew->values()[i] - ans.values()[i]));
            result = 1;
        }
    }

    hvNonsmoothEvalSuccess =
        pnew.evaluateHessianVectorProductNonsmooth(v, quantities);
    ans.setFromFile((char*)"./data/rHv_scaled.txt");
    auto rHvnew = pnew.hessianVectorProductNonsmooth();
    if ((!hvNonsmoothEvalSuccess) || (rHvnew->length() != v->length()))
    {
        result = 1;
    }
    for (int i = 0; i < rHvnew->length(); i++)
    {
        if (fabs(rHvnew->values()[i] - ans.values()[i]) > 1e-7)
        {
            reporter.printf(R_SOLVER, R_BASIC,
                            "Expected Hv[%d]:%8.5f | Actual Hv[%d]:%8.5f | "
                            "Difference: %20.16f\n",
                            ans.values()[i], i, rHvnew->values()[i], i,
                            fabs(rHvnew->values()[i] - ans.values()[i]));
            result = 1;
        }
    }

    // test proximal gradient
    proxComputeSuccess = pnew.computeProximalGradientUpdate(quantities);
    ans.setFromFile((char*)"./data/trueProx_scaled.txt");
    auto proxnew = pnew.proximalGraidentUpdate();
    if ((!proxComputeSuccess) || (proxnew->length() != x->length()))
    {
        result = 1;
    }
    for (int i = 0; i < proxnew->length(); i++)
    {
        if (fabs(proxnew->values()[i] - ans.values()[i]) > 1e-7)
        {
            reporter.printf(R_SOLVER, R_BASIC,
                            "Expected prox[%d]:%8.5f | Actual prox[%d]:%8.5f | "
                            "Difference: %20.16f\n",
                            ans.values()[i], i, proxnew->values()[i], i,
                            fabs(proxnew->values()[i] - ans.values()[i]));
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

}  // end testPointImplementation

#endif /* __TESTPOINT_HPP__ */
