// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#include "FaRSAParameterUpdatePGStepsize.hpp"

#include <cmath>

#include "FaRSADeclarations.hpp"
#include "FaRSADefinitions.hpp"
#include "FaRSALineSearchBacktracking.hpp"

namespace FaRSA
{
// Add options
void ParameterUpdatePGStepsize::addOptions(Options* options, const Reporter* reporter)
{
    options->addStringOption(
        reporter, "PUPG_method", "DecreaseByFraction",
        "Method used to update the stepsize for the proximal gradient update computation.\n "
        "Default     : DecreaseByFraction."
        "Possible values: DecreaseByFraction, LipschitzEstimate, Heuristic");

    options->addDoubleOption(reporter, "PUPG_decrease_fatcor", 0.8, 0.0, 1.0,
                             "Factor used to  decrease the proximal gradient stepsize.\n"
                             "Default     : 0.8.");
    options->addDoubleOption(reporter, "PUPG_increase_fatcor", 1.1, 1.0, FARSA_DOUBLE_INFINITY,
                             "Factor used to  increase the proximal gradient stepsize.\n"
                             "Default     : 1.1.");
    options->addDoubleOption(
        reporter, "PUPG_upper_bound", 1.0, 0.0, FARSA_DOUBLE_INFINITY,
        "The upper bound on the proximal gradient stepsize. (Mainly for GroupFaRSA)\n"
        "Default     : 1.0.");
}  // end addOptions

// Get options
void ParameterUpdatePGStepsize::getOptions(const Options* options, const Reporter* reporter)
{
    options->valueAsString(reporter, "PUPG_method", method_);
    options->valueAsDouble(reporter, "PUPG_decrease_fatcor", decrease_factor_);
    options->valueAsDouble(reporter, "PUPG_increase_fatcor", increase_factor_);
    options->valueAsDouble(reporter, "PUPG_upper_bound", upper_bound_);
}  // end getOptions

// Initialize
void ParameterUpdatePGStepsize::initialize(const Options* options, Quantities* quantities,
                                           const Reporter* reporter)
{
}

void ParameterUpdatePGStepsize::update(const Options* options, Quantities* quantities,
                                       const Reporter* reporter, Strategies* strategies)
{
    auto ptr = std::static_pointer_cast<LineSearchBacktracking>(strategies->lineSearch());
    bool perform_backtrack = ptr->numberOfBacktrack() > 0;
    if (method_.compare("DecreaseByFraction") == 0)
    {
        if (perform_backtrack)
        {
            quantities->setStepsizeProximalGradient(quantities->stepsizeProximalGradient() *
                                                    decrease_factor_);
        }
    }
    else if (method_.compare("Heuristic") == 0)
    {
        if (perform_backtrack)
        {
            quantities->setStepsizeProximalGradient(quantities->stepsizeProximalGradient() *
                                                    decrease_factor_);
        }
        else
        {
            quantities->setStepsizeProximalGradient(quantities->stepsizeProximalGradient() *
                                                    increase_factor_);
        }
    }
    else if (method_.compare("LipschitzEstimate") == 0)
    {
        /*
        update proxStep size according to the paper
            Curtis-Robinson2019_Article_ExploitingNegativeCurvatureInD
        using equation (19) and the one follows.
        */
        double actual_decrease = quantities->trialIterate()->evaluateObjectiveAll(*quantities) -
                                 quantities->currentIterate()->evaluateObjectiveAll(*quantities);
        auto step_taken = quantities->trialIterate()->vector()->makeNewLinearCombination(
            1.0, -1.0, *quantities->currentIterate()->vector());
        double step_taken_norm_square = step_taken->innerProduct(*step_taken);
        double Lipschitz_estimate_old = 1 / quantities->stepsizeProximalGradient();
        double directional_derivative =
            quantities->currentIterate()->gradientSmooth()->innerProduct(*step_taken);
        double model_decrease =
            directional_derivative + 0.5 * step_taken_norm_square * Lipschitz_estimate_old;
        double temp = Lipschitz_estimate_old +
                      2.0 * (actual_decrease - model_decrease) / step_taken_norm_square;
        // safeguard
        if ((actual_decrease - model_decrease) > 0.0)
        {
            Lipschitz_estimate_old =
                fmax(2.0 * Lipschitz_estimate_old, fmin(1e3 * Lipschitz_estimate_old, temp));
        }
        double Lipschitz_estimate_new = fmax(fmax(1e-3, 1e-3 * Lipschitz_estimate_old), temp);
        quantities->setStepsizeProximalGradient(fmin(1 / Lipschitz_estimate_new, upper_bound_));
    }
}
}  // namespace FaRSA
