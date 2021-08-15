// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#include "FaRSAPoint.hpp"

#include <cmath>
#include <iostream>
namespace FaRSA
{
// Constructor, copy elements from input vector
Point::Point(const std::shared_ptr<FunctionSmooth>    function_smooth,
             const std::shared_ptr<FunctionNonsmooth> function_nonsmooth, std::shared_ptr<Vector> vector, double scale)
    : objective_smooth_evaluated_(false),
      objective_nonsmooth_evaluated_(false),
      gradient_smooth_evaluated_(false),
      gradient_nonsmooth_evaluated_(false),
      proximal_gradient_update_evaluated_(false),
      hessian_vector_product_smooth_evaluated_(false),
      hessian_vector_product_nonsmooth_evaluated_(false),
      scale_(scale),
      norm_gradient_all_indices_working_(FARSA_DOUBLE_INFINITY),
      function_smooth_(function_smooth),
      function_nonsmooth_(function_nonsmooth)
{
    objective_smooth_ = FARSA_DOUBLE_INFINITY;
    objective_nonsmooth_ = FARSA_DOUBLE_INFINITY;
    number_of_variables_ = function_smooth_->numberOfVariables();
    number_of_groups_ = function_nonsmooth_->numberOfGroups();
    // Declare new vector
    std::shared_ptr<Vector> new_vector(new Vector(vector->length()));

    // Set point's vector
    vector_ = new_vector;

    // Copy values to point's vector
    vector_->copy(*vector);

    // Set gradient pointer to null
    gradient_smooth_.reset();
    gradient_nonsmooth_.reset();
    proximal_gradient_update_.reset();
    hessian_vector_product_smooth_.reset();
    hessian_vector_product_nonsmooth_.reset();
    per_group_2norm_.reset();
    per_group_gradient_all_2norm_.reset();
    per_group_proximal_gradient_update_2norm_.reset();
}  // end constructor

// Print
void Point::print(const Reporter* reporter, std::string name) const { vector_->print(reporter, name); }

// Make new Point by adding "scalar1" times this Point's vector to "scalar2"
// times other Vector
std::shared_ptr<Point> Point::makeNewLinearCombination(double scalar1, double scalar2, const Vector& other_vector) const
{
    // Create new Vector
    std::shared_ptr<Vector> new_vector = vector_->makeNewLinearCombination(scalar1, scalar2, other_vector);

    // Create new Point
    std::shared_ptr<Point> new_point(new Point(function_smooth_, function_nonsmooth_, new_vector, scale_));

    // Return
    return new_point;

}  // end makeNewLinearCombination

// Determine scale
void Point::determineScale(Quantities& quantities)
{
    // Assert gradient has been evaluated
    ASSERT_EXCEPTION(gradient_smooth_evaluated_, FARSA_GRADIENT_EVALUATION_ASSERT_EXCEPTION,
                     "Gradient should have been evaluated, but wasn't.");

    // Set scale
    if (gradient_smooth_->normInf() > quantities.scalingThreshold())
    {
        scale_ = quantities.scalingThreshold() / gradient_smooth_->normInf();
    }
    quantities.setScaleApplied(scale_);
    // re-scale the related values evaluated point x, where gradient_smooth is
    // evaluated at, since scale_ is the actual one being used instead of 1.0
    if (fabs(1.0 - scale_) > 1e-16)
    {
        scaleObjectiveSmooth();
        scaleGradientSmooth();
        function_nonsmooth_->scalePenalty(scale_);
    }

}  // end determineScale

// Evaluate objective Smooth Part
bool Point::evaluateObjectiveSmooth(Quantities& quantities)
{
    // Check if objective has been evaluated already
    if (!objective_smooth_evaluated_)
    {
        // Set evaluation start time as current time
        clock_t start_time = clock();

        // Evaluate objective value for smooth function
        // deference a shared pointer
        // https://stackoverflow.com/questions/11347111/dereferencing-a-pointer-when-passing-by-reference
        objective_smooth_evaluated_ = function_smooth_->evaluateObjective(*vector_, objective_smooth_);

        // Increment evaluation time
        quantities.incrementEvaluationTime(clock() - start_time);

        // Scale
        objective_smooth_ = scale_ * objective_smooth_;

        // Check for nan
        if (isnan(objective_smooth_))
        {
            objective_smooth_evaluated_ = false;
        }

        // Increment function evaluation counter
        quantities.incrementFunctionCounter();

        // Check for function evaluation limit
        if (quantities.functionCounter() >= quantities.functionEvaluationLimit())
        {
            THROW_EXCEPTION(FARSA_FUNCTION_EVALUATION_LIMIT_EXCEPTION, "Function evaluation limit reached.");
        }

    }  // end if

    // Return
    return objective_smooth_evaluated_;

}  // end evaluateObjective

// Evaluate objective Nonsmooth Part
bool Point::evaluateObjectiveNonsmooth(Quantities& quantities)
{
    // Check if objective has been evaluated already
    if (!objective_nonsmooth_evaluated_)
    {
        // Set evaluation start time as current time
        clock_t start_time = clock();

        // Evaluate objective value
        // deference a shared pointer
        // https://stackoverflow.com/questions/11347111/dereferencing-a-pointer-when-passing-by-reference
        objective_nonsmooth_evaluated_ = function_nonsmooth_->evaluateObjective(*vector_, objective_nonsmooth_);

        // Increment evaluation time
        quantities.incrementEvaluationTime(clock() - start_time);

        // Scale is already done when call determin the scale
        // this is becuase of the proximal gradient update
        // is not a linear transformation, one has to scale
        // penalty first

        // Check for nan
        if (isnan(objective_nonsmooth_))
        {
            objective_nonsmooth_evaluated_ = false;
        }

        /*
        I intentionally comment this out as I think a smooth evalulation
        plus a non evaluation should be counted as one function evaluation
        I don't think the algorithm will only try evaluate nonsmooth
        function at x without evaluating the smooth funcion at x at the
        same time.
        */
        // Increment function evaluation counter
        // quantities.incrementFunctionCounter();

        // // Check for function evaluation limit
        // if (quantities.functionCounter() >=
        //     quantities.functionEvaluationLimit())
        // {
        //     THROW_EXCEPTION(FARSA_FUNCTION_EVALUATION_LIMIT_EXCEPTION,
        //                     "Function evaluation limit reached.");
        // }

    }  // end if

    // Return
    return objective_nonsmooth_evaluated_;

}  // end evaluateObjective

// evaluateObjectiveAll
bool Point::evaluateObjectiveAll(Quantities& quantities)
{
    bool evaluation_success_smooth = evaluateObjectiveSmooth(quantities);
    bool evaluation_success_nonsmooth = evaluateObjectiveNonsmooth(quantities);
    return (evaluation_success_smooth && evaluation_success_nonsmooth);
}  // end evaluateObjectiveAll

// evaluateGradientSmooth
bool Point::evaluateGradientSmooth(Quantities& quantities)
{
    // Check if gradient has been evaluated already
    if ((!gradient_smooth_evaluated_))
    {
        ASSERT_EXCEPTION(objective_smooth_evaluated_, FARSA_FUNCTION_EVALUATION_ASSERT_EXCEPTION,
                         "(From Point::evaluateGradientSmooth): Function Smooth Objective value "
                         "should have been "
                         "evaluated before evaluating the objective, but wasn't.");
        // Declare gradient vector
        std::shared_ptr<Vector> gradient(new Vector(vector_->length()));

        // Set gradient vector
        gradient_smooth_ = gradient;

        // Declare temporary array
        // double* g = new double[vector_->length()];
        Vector g(vector_->length());

        // Set evaluation start time as current time
        clock_t start_time = clock();

        // Evaluate gradient value at the full space
        std::vector<int> full_indicies_vector(vector_->length());
        gradient_smooth_evaluated_ = function_smooth_->evaluateGradient(*vector_, full_indicies_vector, g);

        // Increment evaluation time
        quantities.incrementEvaluationTime(clock() - start_time);

        // Evaluate gradient value
        gradient_smooth_->copyArray(g.values());

        // Scale
        gradient_smooth_->scale(scale_);

        // Check for nan this can be skipped as the implementation for
        // function_smooth_->evaluateGradient will check NA.
        // for (int i = 0; i < gradient_smooth_->length(); i++)
        // {
        //     if (isnan(gradient_smooth_->values()[i]))
        //     {
        //         gradient_smooth_evaluated_ = false;
        //     }
        // }

        // Increment gradient evaluation counter
        quantities.incrementGradientCounter();

        // Check for gradient evaluation limit
        if (quantities.gradientCounter() >= quantities.gradientEvaluationLimit())
        {
            THROW_EXCEPTION(FARSA_GRADIENT_EVALUATION_LIMIT_EXCEPTION, "Gradient evaluation limit reached.");
        }

    }  // end if

    // Return
    return gradient_smooth_evaluated_;

}  // end evaluateGradientSmooth

// evaluateGradientNonsmooth
bool Point::evaluateGradientNonsmooth(Quantities& quantities)
{
    // Check if gradient has been evaluated already
    if ((!gradient_nonsmooth_evaluated_))
    {
        ASSERT_EXCEPTION(objective_nonsmooth_evaluated_, FARSA_FUNCTION_EVALUATION_ASSERT_EXCEPTION,
                         "(Point::evaluateGradientNonsmooth): Function Nonsmooth Objective value "
                         "should have been evaluated before evaluating the gradient, but wasn't.");
        // Declare gradient vector
        std::shared_ptr<Vector> gradient(new Vector(quantities.indiciesWorking()->size()));
        // Set gradient vector
        gradient_nonsmooth_ = gradient;

        // Declare temporary array
        // double* g = new double[vector_->length()];
        Vector g(quantities.indiciesWorking()->size());

        // Set evaluation start time as current time
        clock_t start_time = clock();

        // Evaluate gradient value at the full space
        gradient_nonsmooth_evaluated_ =
            function_nonsmooth_->evaluateGradient(*vector_, *(quantities.indiciesWorking()), g);

        // Increment evaluation time
        quantities.incrementEvaluationTime(clock() - start_time);

        // Evaluate gradient value
        gradient_nonsmooth_->copyArray(g.values());

        // Scale
        // Scale is already done when call determine the scale
        // this is becuase of the proximal gradient update
        // is not a linear transformation, one has to scale
        // penalty first
        // gradient_nonsmooth_->scale(scale_);

    }  // end if
    // Return
    return gradient_nonsmooth_evaluated_;

}  // end evaluateGradient

bool Point::computeProximalGradientUpdate(Quantities& quantities)
{
    if ((!proximal_gradient_update_evaluated_))
    {
        ASSERT_EXCEPTION(gradient_smooth_evaluated_, FARSA_GRADIENT_EVALUATION_ASSERT_EXCEPTION,
                         "(From Point::computeProximalGradientUpdate): Function Smooth Gradient "
                         "should have been evaluated before computing the proximal gradient "
                         "update, but wasn't.");
        // Declare gradient vector
        std::shared_ptr<Vector> proxgrad(new Vector(vector_->length()));
        std::shared_ptr<Vector> proxgradstep(new Vector(vector_->length()));

        // Set proximal_gradient_update vector
        proximal_gradient_update_ = proxgrad;
        proximal_gradient_step_ = proxgradstep;
        // Declare temporary array
        Vector p(vector_->length());

        // Set evaluation start time as current time
        clock_t start_time = clock();

        // Evaluate gradient value at the full space
        proximal_gradient_update_evaluated_ =
            function_nonsmooth_->computeProximalGradientUpdate(*vector_, *gradient_smooth_, quantities, p);

        // Increment evaluation time
        quantities.incrementEvaluationTime(clock() - start_time);

        // Compute proximal gradient update and step
        proximal_gradient_update_->copyArray(p.values());
        // s = xprox - x
        proximal_gradient_step_->copy(*proximal_gradient_update_);
        proximal_gradient_step_->addScaledVector(-1.0, *vector_);
        // // Increment gradient evaluation counter
        // quantities.incrementGradientCounter();

        // // Check for gradient evaluation limit
        // if (quantities.gradientCounter() >=
        //     quantities.gradientEvaluationLimit())
        // {
        //     THROW_EXCEPTION(FARSA_GRADIENT_EVALUATION_LIMIT_EXCEPTION,
        //                     "Gradient evaluation limit reached.");
        // }

    }  // end if

    // Return
    return proximal_gradient_update_evaluated_;
}

bool Point::evaluateHessianVectorProductSmooth(std::shared_ptr<Vector> v, Quantities& quantities)
{
    ASSERT_EXCEPTION(gradient_smooth_evaluated_, FARSA_GRADIENT_EVALUATION_ASSERT_EXCEPTION,
                     "(From Point::evaluateHessianVectorProductSmooth): Function Smooth Gradient value "
                     "should have been evaluated before evaluating the Hessian Vector product, but "
                     "wasn't.");
    // Declare gradient vector
    std::shared_ptr<Vector> hv(new Vector(quantities.indiciesWorking()->size()));

    // Set gradient vector
    hessian_vector_product_smooth_ = hv;

    // Declare temporary array
    Vector temp(v->length());

    // Set evaluation start time as current time
    clock_t start_time = clock();

    // Evaluate gradient value at the full space
    hessian_vector_product_smooth_evaluated_ =
        function_smooth_->evaluateHessianVectorProduct(*vector_, *(quantities.indiciesWorking()), *v, temp);

    // Increment evaluation time
    quantities.incrementEvaluationTime(clock() - start_time);

    // Evaluate gradient value
    hessian_vector_product_smooth_->copyArray(temp.values());

    // Scale
    hessian_vector_product_smooth_->scale(scale_);

    // Increment Hessian Vector Product counter
    quantities.incrementHessianVectorCounter();

    // Check for gradient evaluation limit
    if (quantities.hessianVectorCounter() >= quantities.hessianVectorProductEvaluationLimit())
    {
        THROW_EXCEPTION(FARSA_HESSIAN_VECTOR_PRODUCT_EVALUATION_LIMIT_EXCEPTION,
                        "Hessian Vector Product evaluation limit reached.");
    }

    // Return
    return hessian_vector_product_smooth_evaluated_;
}
bool Point::evaluateHessianVectorProductNonsmooth(std::shared_ptr<Vector> v, Quantities& quantities)
{
    ASSERT_EXCEPTION(gradient_nonsmooth_evaluated_, FARSA_GRADIENT_EVALUATION_ASSERT_EXCEPTION,
                     "(From Point::evaluateHessianVectorProductNonsmooth): Function Smooth "
                     "Gradient value should have been "
                     "evaluated before evaluating the Hessian Vector "
                     "product, but wasn't.");
    // Declare gradient vector
    std::shared_ptr<Vector> hv(new Vector(quantities.indiciesWorking()->size()));

    // Set gradient vector
    hessian_vector_product_nonsmooth_ = hv;

    // Declare temporary array
    Vector temp(v->length());

    // Set evaluation start time as current time
    clock_t start_time = clock();

    // Evaluate gradient value at the full space
    hessian_vector_product_nonsmooth_evaluated_ =
        function_nonsmooth_->evaluateHessianVectorProduct(*vector_, *(quantities.indiciesWorking()), *v, temp);

    // Increment evaluation time
    quantities.incrementEvaluationTime(clock() - start_time);

    // Evaluate gradient value
    hessian_vector_product_nonsmooth_->copyArray(temp.values());

    // Scale
    // Scale is already done when call determine the scale
    // this is becuase of the proximal gradient update
    // is not a linear transformation, one has to scale
    // penalty first
    // hessian_vector_product_nonsmooth_->scale(scale_);

    // // Check for gradient evaluation limit
    // if (quantities.hessianVectorCounter() >=
    //     quantities.hessianVectorProductEvaluationLimit())
    // {
    //     THROW_EXCEPTION(
    //         FARSA_HESSIAN_VECTOR_PRODUCT_EVALUATION_LIMIT_EXCEPTION,
    //         "Hessian Vector Product evaluation limit reached.");
    // }

    // Return
    return hessian_vector_product_nonsmooth_evaluated_;
}

bool Point::evaluatePerGroupStatistics(Quantities& quantities)
{
    if ((!per_group_2norm_evaluated_) or (!per_group_gradient_all_2norm_evaluated_) or
        (!per_group_proximal_gradient_update_2norm_evaluated_))
    {
        std::shared_ptr<Vector> per_group_2norm(new Vector(quantities.numberOfGroups()));
        per_group_2norm_ = per_group_2norm;
        std::shared_ptr<Vector> per_group_gradient_all_2norm(new Vector(quantities.numberOfGroups()));
        per_group_gradient_all_2norm_ = per_group_gradient_all_2norm;
        std::shared_ptr<Vector> per_group_proximal_gradient_update_2norm(new Vector(quantities.numberOfGroups()));
        per_group_proximal_gradient_update_2norm_ = per_group_proximal_gradient_update_2norm;
        if (quantities.groupsFormat().compare("C") == 0)
        {
            auto groups = quantities.groups();
            // create working indicies (only need to evaluate gradient over these indices)
            std::shared_ptr<std::vector<int>> indicies_working(new std::vector<int>());

            // compute per_group_2norm_ and per_group_proximal_gradient_update_2norm_ and collect non-zero group
            // indicies
            for (int i = 0; i < quantities.numberOfGroups(); i++)
            {
                double group_i_norm = 0.0;
                double group_i_proximal_gradient_update_norm = 0.0;
                int    start = (*groups)[i][0];
                int    end = (*groups)[i][1];
                for (int j = start; j <= end; j++)
                {
                    group_i_norm += pow(vector_->values()[j], 2);
                    group_i_proximal_gradient_update_norm += pow(proximalGraidentUpdate()->values()[j], 2);
                }
                per_group_2norm_->valuesModifiable()[i] = sqrt(group_i_norm);
                per_group_proximal_gradient_update_2norm_->valuesModifiable()[i] =
                    sqrt(group_i_proximal_gradient_update_norm);
                per_group_gradient_all_2norm_->valuesModifiable()[i] = FARSA_DOUBLE_INFINITY;
                // current group is non-zero group add all indices to the indicies vector
                if (fabs(group_i_norm - 0.0) > 1e-16)
                {
                    for (int j = start; j <= end; j++)
                    {
                        indicies_working->push_back(j);
                    }
                }
            }

            // compute per_group_gradient_all_2norm_
            quantities.setIndiciesWorking(indicies_working);
            evaluateGradientSmooth(quantities);
            evaluateGradientNonsmooth(quantities);
            for (int i = 0; i < quantities.numberOfGroups(); i++)
            {
                // current group is non-zero group
                if (fabs(per_group_2norm_->values()[i] - 0.0) > 1e-16)
                {
                    int    start = (*groups)[i][0];
                    int    end = (*groups)[i][1];
                    double group_i_gradient_all_norm = 0.0;
                    for (int j = start; j <= end; j++)
                    {
                        group_i_gradient_all_norm +=
                            pow(gradient_smooth_->values()[j] + gradient_nonsmooth_->values()[j], 2);
                    }
                    per_group_gradient_all_2norm_->valuesModifiable()[i] = sqrt(group_i_gradient_all_norm);
                }
            }
            per_group_2norm_evaluated_ = true;
            per_group_gradient_all_2norm_evaluated_ = true;
            per_group_proximal_gradient_update_2norm_evaluated_ = true;
        }
        else if (quantities.groupsFormat().compare("NC") == 0)
        {
            THROW_EXCEPTION(FARSA_NO_IMPLEMENTATION_EXCEPTION, "No implementation for the group_format_!");
        }
        else
        {
            THROW_EXCEPTION(FARSA_NO_IMPLEMENTATION_EXCEPTION, "No implementation for the group_format_!");
        }
    }
    // Return
    return per_group_2norm_evaluated_ && per_group_gradient_all_2norm_evaluated_ &&
           per_group_proximal_gradient_update_2norm_evaluated_;
}
}  // namespace FaRSA
