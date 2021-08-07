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
    GroupL1(std::shared_ptr<std::vector<std::vector<int>>> groups,
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
    double const                         penalty() const { return penalty_; };
    std::shared_ptr<std::vector<double>> weights() { return weights_; };
    std::shared_ptr<std::vector<std::vector<int>>> const groups() const
    {
        return groups_;
    };
    //@}

    /** @name Evaluate methods */
    //@{
    bool scalePenalty(double scale);
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
    bool evaluateGradient(const Vector& x, const std::vector<int>& indicies,
                          Vector& g);
    /**
     * Evaluates gradient in subspace defined by variables
     * in subgroups.
     * \param[in] x is a given point/iterate, a constant double array
     * \param[in] indicies is a vector of column indices
     * \param[in] v is a given vector, a constant double array
     * \param[out] Hv is the product of the Hessian and "v", a double array
     * (return value) \return indicator of success (true) or failure (false)
     */
    bool evaluateHessianVectorProduct(const Vector&           x,
                                      const std::vector<int>& indicies,
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
    std::shared_ptr<std::vector<double>> weights_;
    /** a vector of vector specifing the groups of variables
     *  For example, consider the following group structures:
     *  group1: [0,1]; group2: [2,3,4,5]; group3:[6]
     *  the groups_ [[0,1], [2,5], [6,6]], for i-th subvector,
     *  the first element is the statring position of the i-th group
     *  the second element is the ending position of the i-th group
     */
    std::vector<double> per_group_norm_;
    // std::vector<int>              per_group_size_;
    std::vector<int>      idx_to_group_;
    std::map<int, double> per_group_gradient_step_norm_;
    //@}

};  // end GroupL1

#endif /* __GROUPL1_HPP__ */
