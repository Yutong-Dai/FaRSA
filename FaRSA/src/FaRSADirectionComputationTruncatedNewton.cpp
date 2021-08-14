// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#include "FaRSADirectionComputationTruncatedNewton.hpp"

#include "FaRSADeclarations.hpp"
#include "FaRSADefinitions.hpp"

namespace FaRSA
{
// Add options
void DirectionComputationTruncatedNewton::addOptions(Options* options, const Reporter* reporter)
{
    options->addDoubleOption(reporter, "DCTN_maxCGIters", FARSA_DOUBLE_INFINITY, 0.0, FARSA_DOUBLE_INFINITY,
                             "Max number of conjugate gradient steps allowed to perform.\n"
                             "Default     : FARSA_DOUBLE_INFINITY.");
    options->addDoubleOption(reporter, "DCTN_CGBigFactor", 1e3, 0.0, FARSA_DOUBLE_INFINITY,
                             "Max number of conjugate gradient steps allowed to perform.\n"
                             "Default     : 1e3.");
    options->addBoolOption(reporter, "DCTN_verbose", true,
                           "A parameter controls whether should print more details.\n"
                           "Default     : true.");
}  // end addOptions

// Set options
void DirectionComputationTruncatedNewton::getOptions(const Options* options, const Reporter* reporter)
{
    options->valueAsDouble(reporter, "DCTN_maxCGIters", max_CG_iters_);
    options->valueAsDouble(reporter, "DCTN_CGBigFactor", cg_big_factor_);
    options->valueAsBool(reporter, "DCTN_verbose", verbose_);
}  // end getOptions

// Initialize
void DirectionComputationTruncatedNewton::initialize(const Options* options, Quantities* quantities,
                                                     const Reporter* reporter)
{
}

// Compute direction
void DirectionComputationTruncatedNewton::computeDirection(const Options* options, Quantities* quantities,
                                                           const Reporter* reporter, Strategies* strategies)
{
    // Initialize values
    setStatus(DC_UNSET);
    // if perform computation
    if (performCompuation())
    {
        // try direction computation, terminate on any exception
        try
        {
            // Initialize boolean for evaluation
            bool evaluation_success = false;

            // Evaluate current gradient smooth
            evaluation_success = quantities->currentIterate()->evaluateGradientSmooth(*quantities);

            // Check for successful evaluation
            if (!evaluation_success)
            {
                THROW_EXCEPTION(DC_EVALUATION_FAILURE_EXCEPTION,
                                "Direction computation unsuccessful. "
                                "Gradient smooth evaluation failed.");
            }
            // Evaluate current gradient nonsmooth
            evaluation_success = quantities->currentIterate()->evaluateGradientNonsmooth(*quantities);

            // Check for successful evaluation
            if (!evaluation_success)
            {
                THROW_EXCEPTION(DC_EVALUATION_FAILURE_EXCEPTION,
                                "Direction computation unsuccessful. "
                                "Gradient nonsmoth evaluation failed.");
            }
            // form the gradient of smooth and nonsmooth at the subspace
            int    dimension = quantities->indiciesWorking()->size();
            Vector gradF_indices_working(dimension);
            auto   current_iterate = quantities->currentIterate();
            for (int i = 0; i < dimension; i++)
            {
                auto idx = (*quantities->indiciesWorking())[i];
                gradF_indices_working.valuesModifiable()[i] = current_iterate->gradientSmooth()->values()[idx] +
                                                              current_iterate->gradientNonsmooth()->values()[idx];
            }
            // initial residual
            Vector residual(dimension);
            residual.copy(gradF_indices_working);
            double norm_residual = sqrt(residual.innerProduct(residual));

            // initialize conjugate direction p as the negative gradient
            std::shared_ptr<Vector> p(new Vector(dimension));
            p->addScaledVector(-1.0, gradF_indices_working);

            // starting point
            Vector d(dimension);

            // set up termination conditions
            int    max_iters = fmin(dimension, max_CG_iters_);
            double target_residual = fmax(1e-10, fmin(0.1, pow(norm_residual, 0.5)) * norm_residual);
            double norm_gradF_indices_working = sqrt(gradF_indices_working.innerProduct(gradF_indices_working));

            // CG iterations
            int                     iterations = 0;
            double                  pTHp = 0.0;
            double                  alphaCG = 0.0;
            double                  betaCG = 0.0;
            double                  norm_d = 0.0;
            double                  norm_residual_old = 0.0;
            std::shared_ptr<Vector> Hp;
            std::string             flag;
            bool                    termination = false;
            while (!termination)
            {
                iterations += 1;
                evaluation_success = current_iterate->evaluateHessianVectorProductSmooth(p, *quantities);
                if (!evaluation_success)
                {
                    THROW_EXCEPTION(DC_EVALUATION_FAILURE_EXCEPTION,
                                    "Direction computation unsuccessful. "
                                    "Hessian vector product smooth evaluation failed.");
                }
                evaluation_success = current_iterate->evaluateHessianVectorProductNonsmooth(p, *quantities);
                if (!evaluation_success)
                {
                    THROW_EXCEPTION(DC_EVALUATION_FAILURE_EXCEPTION,
                                    "Direction computation unsuccessful. "
                                    "Hessian vector product nonsmooth evaluation failed.");
                }
                Hp = current_iterate->hessianVectorProductSmooth()->makeNewCopy();
                Hp->addScaledVector(1.0, *current_iterate->hessianVectorProductNonsmooth());
                pTHp = p->innerProduct(*Hp);
                alphaCG = pow(norm_residual, 2.0) / pTHp;
                d.addScaledVector(alphaCG, *p);
                norm_d = sqrt(d.innerProduct(d));
                residual.addScaledVector(alphaCG, *Hp);
                norm_residual_old = norm_residual;
                // update norm residual
                norm_residual = sqrt(residual.innerProduct(residual));

                // check termination
                if (norm_residual <= target_residual)
                {
                    flag = "CGtol";
                    termination = true;
                }
                else if (norm_d > cg_big_factor_ * fmin(1, norm_gradF_indices_working))
                {
                    flag = "CGbig";
                    termination = true;
                }
                else if (iterations > max_iters)
                {
                    flag = "CGmax";
                    termination = true;
                }
                // compute next conjugate direction
                betaCG = pow(norm_residual, 2.0) / pow(norm_residual_old, 2.0);
                p->linearCombination(-1.0, residual, betaCG, *p);
            }

            Vector search_direction_actual(quantities->numberOfVariables());
            // filling the second order direction
            for (int i = 0; i < dimension; i++)
            {
                int idx = (*quantities->indiciesWorking())[i];
                search_direction_actual.valuesModifiable()[idx] = d.values()[i];
            }
            quantities->directionSecondOrder()->copy(search_direction_actual);

            // Set status
            setStatus(DC_SUCCESS);
            // print information
            if (verbose_)
            {
                reporter->printf(R_SOLVER, R_PER_ITERATION, " %6d %+.2e %s %5d %+.2e %+.2e", dimension,
                                 norm_gradF_indices_working, flag.c_str(), iterations, norm_residual, target_residual);
            }
            // Check for success
            THROW_EXCEPTION(DC_SUCCESS_EXCEPTION, "Direction computation successful.")

        }  // end try

        // catch exceptions
        catch (DC_SUCCESS_EXCEPTION& exec)
        {
            setStatus(DC_SUCCESS);
        }
        catch (DC_EVALUATION_FAILURE_EXCEPTION& exec)
        {
            setStatus(DC_EVALUATION_FAILURE);
        }
    }
    else
    {
        setStatus(DC_SKIPPED);
        if (verbose_)
        {
            reporter->printf(R_SOLVER, R_PER_ITERATION, " ------ --------- ----- ----- --------- ---------");
        }
    }

}  // end computeDirection

}  // namespace FaRSA
