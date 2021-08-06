// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

// Description : Implementation for Group-L1
//                 r(x) = penalty * \sum_{i=1}^K w_i||x_{Gi}||_2

#ifndef __GROUPL1_HPP__
#define __GROUPL1_HPP__

#include <cassert>
#include <cmath>
#include <map>
#include <vector>

#include "FaRSAFunctionNonsmooth.hpp"
using namespace FaRSA;

class GroupL1 : public FunctionNonsmooth
{
   public:
    /** @name Constructors */
    //@{
    /**
     * Constructor
     */
    GroupL1(std::shared_ptr<std::vector<std::vector<int>>> group,
            std::shared_ptr<std::vector<double>> weights, double penalty);
    //@}

    /** @name Destructor */
    //@{
    /**
     * Destructor
     */
    ~GroupL1(){};

    /** @name Get methods */
    //@{
    double                               penalty() const { return penalty_; };
    std::shared_ptr<std::vector<double>> weights() const { return weights_; };
    //@}

    /** @name Evaluate methods */
    //@{
    /**
     * Evaluates objective
     * \param[in] x is a given point/iterate, a constant double array
     * \param[out] f is the objective value at "x", a double (return value)
     * \return indicator of success (true) or failure (false)
     */
    bool evaluateObjective(const Vector& x, double& f);
    /**
     * Evaluates gradient
     * \param[in] x is a given point/iterate, a constant double array
     * \param[out] g is the gradient value at "x", a double array (return value)
     * \return indicator of success (true) or failure (false)
     */
    bool evaluateGradient(const Vector& x, const std::vector<int>& group_idx,
                          Vector& g);
    /**
     * Evaluates gradient in subspace defined by variables
     * in subgroups.
     * \param[in] x is a given point/iterate, a constant double array
     * \param[in] group_idx is a vector of column indices
     * \param[in] v is a given vector, a constant double array
     * \param[out] Hv is the product of the Hessian and "v", a double array
     * (return value) \return indicator of success (true) or failure (false)
     */
    bool evaluateHessianVectorProduct(const Vector&           x,
                                      const std::vector<int>& group_idx,
                                      const Vector& v, Vector& Hv);
    /**
     * @brief Evaluator proximal operator at u
     *
     * \param u
     * \param quantities
     * \return true
     * \return false
     */
    bool computeProximalGradientUpdate(const Vector& x, const Vector& gradfx,
                                       Quantities& quantities,
                                       Vector&     proxgrad);

    //@}

   private:
    /** @name Default compiler generated methods
     * (Hidden to avoid implicit creation/calling.)
     */
    //@{
    /**
     * Copy constructor
     */
    GroupL1(const GroupL1&);
    /**
     * Overloaded equals operator
     */
    void operator=(const GroupL1&);
    //@}
    /**
     * Private members
     */
    //@{
    double                               penalty_;
    std::shared_ptr<std::vector<double>> weights_;

    std::vector<double> per_group_norm_;
    // std::vector<int>              per_group_size_;
    std::vector<int>      idx_to_group_;
    std::map<int, double> per_group_gradient_step_norm_;
    //@}

};  // end GroupL1

#endif /* __GROUPL1_HPP__ */
