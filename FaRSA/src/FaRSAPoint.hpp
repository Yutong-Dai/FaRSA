// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __FARSAPOINT_HPP__
#define __FARSAPOINT_HPP__

#include <memory>
#include <string>

#include "FaRSADeclarations.hpp"
#include "FaRSAFunctionNonsmooth.hpp"
#include "FaRSAFunctionSmooth.hpp"
#include "FaRSAQuantities.hpp"
#include "FaRSAReporter.hpp"
#include "FaRSAVector.hpp"
namespace FaRSA
{
/**
 * Forward declarations
 */
class Quantities;
class FunctionSmooth;
class FunctionNonsmooth;
class Reporter;
class Vector;

/**
 * Point class
 */
class Point
{
   public:
    /** @name Constructors */
    //@{
    /**
     * Construct Point
     * \param[in] f is pointer to FaRSAFunctionSmooth object
     * \param[in] r is pointer to FaRSAFunctionNonsmooth object
     * \param[in] vector is pointer to Vector object
     * \param[in] scale is scaling factor for function evaluations
     */
    Point(const std::shared_ptr<FunctionSmooth> f, const std::shared_ptr<FunctionNonsmooth> r,
          const std::shared_ptr<Vector> vector, double scale);
    //@}

    /** @name Destructor */
    //@{
    /**
     * Delete data
     */
    ~Point(){};
    //@}

    /** @name Print methods */
    //@{
    /**
     * Print data
     * \param[in] reporter is pointer to Reporter object from FaRSA
     * \param[in] name is name of point to print
     */
    void print(const Reporter* reporter, std::string name) const;
    //@}

    /** @name Make-new methods */
    //@{
    /**
     * Make new Point by adding "scalar1" times this Point's Vector to "scalar2"
     * times other_vector \param[in] scalar1 is scalar value for linear
     * combination \param[in] scalar2 is scalar value for linear combination
     * \param[in] other_vector is reference to other Vector
     * \return pointer to new Point
     */
    std::shared_ptr<Point> makeNewLinearCombination(double scalar1, double scalar2,
                                                    const Vector& other_vector) const;
    //@}

    /** @name Set methods */
    //@{
    /**
     * Determine scale for objective function
     * \param[in,out] quantities is reference to IterationQuanitites object from
     * FaRSA \return boolean indicating success
     */
    void determineScale(Quantities& quantities);
    /**
     * Evaluate objective
     * \param[in,out] quantities is reference to IterationQuanitites object from
     * FaRSA \return boolean indicating success
     */
    bool evaluateObjectiveSmooth(Quantities& quantities);
    /**
     * Evaluate objective
     * \param[in,out] quantities is reference to IterationQuanitites object from
     * FaRSA \return boolean indicating success
     */
    bool evaluateObjectiveNonsmooth(Quantities& quantities);
    /**
     * Evaluate gradient
     * \param[in,out] quantities is reference to IterationQuanitites object from
     * FaRSA \return boolean indicating success
     */
    bool evaluateObjectiveAll(Quantities& quantities);
    bool evaluateGradientSmooth(Quantities& quantities);
    /**
     * Evaluate gradient
     * \param[in,out] quantities is reference to IterationQuanitites object from
     * FaRSA \return boolean indicating success
     */
    bool evaluateGradientNonsmooth(Quantities& quantities);
    bool evaluateHessianVectorProductSmooth(std::shared_ptr<Vector> v, Quantities& quantities);
    bool evaluateHessianVectorProductNonsmooth(std::shared_ptr<Vector> v, Quantities& quantities);

    // set both proximal gradient update and proximal gradient step
    bool computeProximalGradientUpdate(Quantities& quantities);

    /**
     * Scale objective
     */
    void scaleObjectiveSmooth()
    {
        ASSERT_EXCEPTION(objective_smooth_evaluated_, FARSA_FUNCTION_EVALUATION_ASSERT_EXCEPTION,
                         "Objective Smooth should have been evaluated, but wasn't.");
        objective_smooth_ *= scale_;
    }
    void scaleObjectiveNonsmooth()
    {
        ASSERT_EXCEPTION(objective_nonsmooth_evaluated_, FARSA_FUNCTION_EVALUATION_ASSERT_EXCEPTION,
                         "Objective Nonsmooth should have been evaluated, but wasn't.");
        objective_nonsmooth_ *= scale_;
    }
    /**
     * Scale gradient
     */
    void scaleGradientSmooth()
    {
        ASSERT_EXCEPTION(gradient_smooth_evaluated_, FARSA_GRADIENT_EVALUATION_ASSERT_EXCEPTION,
                         "Gradient Smooth should have been evaluated, but wasn't.");
        gradient_smooth_->scale(scale_);
    };
    void scaleGradientNonsmooth()
    {
        ASSERT_EXCEPTION(gradient_nonsmooth_evaluated_, FARSA_GRADIENT_EVALUATION_ASSERT_EXCEPTION,
                         "Gradient Nonsmooth should have been evaluated, but wasn't.");
        gradient_nonsmooth_->scale(scale_);
    };
    //@}

    /** @name Get methods */
    //@{
    /**
     * Get pointer to Problem
     * \return is pointer to Problem corresponding to point
     */
    inline std::shared_ptr<FunctionSmooth>    funcion_smooth() const { return function_smooth_; };
    inline std::shared_ptr<FunctionNonsmooth> funcion_nonsmooth() const
    {
        return function_nonsmooth_;
    };
    /**
     * Get pointer to vector
     * \return is pointer to Vector defining point
     */
    inline std::shared_ptr<Vector> vector() const { return vector_; };
    /**
     * Get objective
     * \return objective as double
     */
    inline double objectiveSmooth() const
    {
        ASSERT_EXCEPTION(objective_smooth_evaluated_, FARSA_FUNCTION_EVALUATION_ASSERT_EXCEPTION,
                         "FunctionSmooth should have been evaluated, but wasn't.");
        return objective_smooth_;
    };

    inline double objectiveNonsmooth() const
    {
        ASSERT_EXCEPTION(objective_nonsmooth_evaluated_, FARSA_FUNCTION_EVALUATION_ASSERT_EXCEPTION,
                         "FunctionNonsmooth should have been evaluated, but wasn't.");
        return objective_nonsmooth_;
    };

    inline double objectiveAll() const
    {
        ASSERT_EXCEPTION(objective_smooth_evaluated_, FARSA_FUNCTION_EVALUATION_ASSERT_EXCEPTION,
                         "FunctionSmooth should have been evaluated, but wasn't.");
        ASSERT_EXCEPTION(objective_nonsmooth_evaluated_, FARSA_FUNCTION_EVALUATION_ASSERT_EXCEPTION,
                         "FunctionNonsmooth should have been evaluated, but wasn't.");
        return objective_smooth_ + objective_nonsmooth_;
    }

    /**
     * Get objective (unscaled)
     * \return objective (unscaled) as double
     */
    inline double objectiveSmoothUnscaled() const
    {
        ASSERT_EXCEPTION(objective_smooth_evaluated_, FARSA_FUNCTION_EVALUATION_ASSERT_EXCEPTION,
                         "Function Smooth should have been evaluated, but wasn't.");
        return objective_smooth_ / scale_;
    };

    inline double objectiveNonsmoothUnscaled() const
    {
        ASSERT_EXCEPTION(objective_nonsmooth_evaluated_, FARSA_FUNCTION_EVALUATION_ASSERT_EXCEPTION,
                         "Function Nonsmooth should have been evaluated, but wasn't.");
        return objective_nonsmooth_ / scale_;
    };
    /**
     * Get pointer to gradient
     * \return gradient as pointer to Vector
     */
    inline std::shared_ptr<Vector> gradientSmooth() const
    {
        ASSERT_EXCEPTION(gradient_smooth_evaluated_, FARSA_GRADIENT_EVALUATION_ASSERT_EXCEPTION,
                         "Gradient Smooth should have been evaluated, but wasn't.");
        return gradient_smooth_;
    };
    inline std::shared_ptr<Vector> gradientNonsmooth() const
    {
        ASSERT_EXCEPTION(gradient_smooth_evaluated_, FARSA_GRADIENT_EVALUATION_ASSERT_EXCEPTION,
                         "Gradient Nonsmooth should have been evaluated, but wasn't.");
        return gradient_nonsmooth_;
    };
    inline std::shared_ptr<Vector> proximalGraidentUpdate() const
    {
        ASSERT_EXCEPTION(proximal_gradient_update_evaluated_,
                         FARSA_PROXIMAL_GRADIENT_COMPUTATION_ASSERT_EXCEPTION,
                         "Proximal Gradient Update should have been "
                         "computed, but wasn't.");
        return proximal_gradient_update_;
    };
    inline std::shared_ptr<Vector> proximalGraidentStep() const
    {
        ASSERT_EXCEPTION(proximal_gradient_update_evaluated_,
                         FARSA_PROXIMAL_GRADIENT_COMPUTATION_ASSERT_EXCEPTION,
                         "Proximal Gradient Update should have been "
                         "computed, but wasn't.");
        return proximal_gradient_step_;
    };

    inline std::shared_ptr<Vector> hessianVectorProductSmooth() const
    {
        ASSERT_EXCEPTION(hessian_vector_product_smooth_evaluated_,
                         FARSA_HESSIAN_VECTOR_PRODUCT_EVALUATION_ASSERT_EXCEPTION,
                         "Hessian Vector Product Smooth part should have been "
                         "evaluated, but wasn't.");
        return hessian_vector_product_smooth_;
    };

    inline std::shared_ptr<Vector> hessianVectorProductNonsmooth() const
    {
        ASSERT_EXCEPTION(hessian_vector_product_nonsmooth_evaluated_,
                         FARSA_HESSIAN_VECTOR_PRODUCT_EVALUATION_ASSERT_EXCEPTION,
                         "Hessian Vector Product Nonsmooth part should have been "
                         "evaluated, but wasn't.");
        return hessian_vector_product_nonsmooth_;
    };
    /**
     * Get scale
     * \return is scale factor
     */
    inline double scale() const { return scale_; };
    //@}

   private:
    /** @name Default compiler generated methods
     * (Hidden to avoid implicit creation/calling.)
     */
    //@{
    /**
     * Copy constructor
     */
    Point(const Point&);
    /**
     * Overloaded equals operator
     */
    void operator=(const Point&);
    //@}

    /** @name Private members */
    //@{
    bool objective_smooth_evaluated_;
    bool objective_nonsmooth_evaluated_;
    bool gradient_smooth_evaluated_;
    bool gradient_nonsmooth_evaluated_;
    bool proximal_gradient_update_evaluated_;
    bool hessian_vector_product_smooth_evaluated_;
    bool hessian_vector_product_nonsmooth_evaluated_;

    double                             objective_smooth_;
    double                             objective_nonsmooth_;
    double                             scale_;
    int                                number_of_variables_;
    int                                number_of_groups_;
    std::shared_ptr<FunctionSmooth>    function_smooth_;
    std::shared_ptr<FunctionNonsmooth> function_nonsmooth_;
    std::shared_ptr<Vector>            vector_;
    std::shared_ptr<Vector>            gradient_smooth_;
    std::shared_ptr<Vector>            gradient_nonsmooth_;
    std::shared_ptr<Vector>            proximal_gradient_update_;
    std::shared_ptr<Vector>            proximal_gradient_step_;
    std::shared_ptr<Vector>            hessian_vector_product_smooth_;
    std::shared_ptr<Vector>            hessian_vector_product_nonsmooth_;
    //@}

};  // end Point

}  // namespace FaRSA

#endif /* __FARSAPOINT_HPP__ */
