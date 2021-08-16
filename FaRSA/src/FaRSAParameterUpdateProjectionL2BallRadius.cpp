// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#include "FaRSAParameterUpdateProjectionL2BallRadius.hpp"

#include <cmath>
#include <iostream>

#include "FaRSADeclarations.hpp"
#include "FaRSADefinitions.hpp"
#include "FaRSALineSearchBacktracking.hpp"
#include "FaRSASpacePartitionGroupL1PGBased.hpp"

namespace FaRSA
{
// Add options
void ParameterUpdateProjectionL2BallRadius::addOptions(Options* options, const Reporter* reporter)
{
    options->addDoubleOption(reporter, "PU_kappa1_max", 1e+06, 0.0, FARSA_DOUBLE_INFINITY,
                             "kappa1_max.\n"
                             "Default     : 1e+06.");
    options->addDoubleOption(reporter, "PU_kappa1_min", 1e-05, 0.0, FARSA_DOUBLE_INFINITY,
                             "kappa1_min.\n"
                             "Default     : 1e-05.");
    options->addDoubleOption(reporter, "PU_kappa2_max", 1e+05, 0.0, FARSA_DOUBLE_INFINITY,
                             "kappa2_max.\n"
                             "Default     : 1e+05.");
    options->addDoubleOption(reporter, "PU_kappa2_min", 1e-06, 0.0, FARSA_DOUBLE_INFINITY,
                             "kappa2_min.\n"
                             "Default     : 1e-06.");
    options->addDoubleOption(reporter, "PU_kappa_increase_factor", 10.0, 1.0, FARSA_DOUBLE_INFINITY,
                             "Factor to increase the kappa1 and kappa2.\n"
                             "Default     : 10.");
    options->addDoubleOption(reporter, "PU_kappa_decrease_factor", 0.1, 0.0, 1.0,
                             "Factor to decrease the kappa1 and kappa2.\n"
                             "Default     : 0.1.");
}  // end addOptions

// Get options
void ParameterUpdateProjectionL2BallRadius::getOptions(const Options* options, const Reporter* reporter)
{
    options->valueAsDouble(reporter, "PU_kappa1_max", kappa1_max_);
    options->valueAsDouble(reporter, "PU_kappa1_min", kappa1_min_);
    options->valueAsDouble(reporter, "PU_kappa2_max", kappa2_max_);
    options->valueAsDouble(reporter, "PU_kappa2_min", kappa2_min_);
    options->valueAsDouble(reporter, "PU_kappa_increase_factor", kappa_increase_factor_);
    options->valueAsDouble(reporter, "PU_kappa_decrease_factor", kappa_decrease_factor_);
}  // end getOptions

// Initialize
void ParameterUpdateProjectionL2BallRadius::initialize(const Options* options, Quantities* quantities,
                                                       const Reporter* reporter)
{
}

void ParameterUpdateProjectionL2BallRadius::update(const Options* options, Quantities* quantities,
                                                   const Reporter* reporter, Strategies* strategies)
{
    // Initialize values
    setStatus(PU_UNSET);
    auto ptr_sp = std::static_pointer_cast<SpacePartitionGroupL1PGBased>(strategies->spacePartition());
    auto ptr_ls = std::static_pointer_cast<LineSearchBacktracking>(strategies->lineSearch());
    try
    {
        if (ptr_ls->iterationType().compare("first_order") == 0)
        {
            if (quantities->consequtiveFirstOrderIterationCounter() > 5)
            {
                ptr_sp->setKappa1(fmax(ptr_sp->kappa1() * kappa_decrease_factor_, kappa1_min_));
                ptr_sp->setKappa2(fmax(ptr_sp->kappa2() * kappa_decrease_factor_, kappa2_min_));
            }
        }
        else
        {
            if (ptr_ls->numberOfBacktrack() + ptr_ls->numberOfProjection() > 6)
            {
                ptr_sp->setKappa1(fmin(ptr_sp->kappa1() * kappa_increase_factor_, kappa1_max_));
                ptr_sp->setKappa2(fmin(ptr_sp->kappa2() * kappa_increase_factor_, kappa2_max_));
            }
        }
    }
    catch (PU_SUCCESS_EXCEPTION& exec)
    {
        setStatus(PU_SUCCESS);
    }
    catch (PU_EVALUATION_FAILURE_EXCEPTION& exec)
    {
        setStatus(PU_EVALUATION_FAILURE);
    }
}
}  // namespace FaRSA
