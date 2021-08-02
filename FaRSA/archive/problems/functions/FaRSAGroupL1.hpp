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

#include "FaRSAFunctionNonSmooth.hpp"
using namespace FaRSA;

class GroupL1 : public FunctionNonSmooth
{
   public:
    /** @name Constructors */
    //@{
    /**
     * Constructor
     */
    GroupL1(std::vector<std::vector<int>>& groups, std::vector<float>& weights,
            float penalty);
    //@}

    /** @name Destructor */
    //@{
    /**
     * Destructor
     */
    ~GroupL1(){};

    /** @name Get methods */
    //@{
    float                          getPenalty() const { return penalty_; };
    std::vector<float>&            getWeights() { return weights_; };
    std::vector<std::vector<int>>& getGroups() { return groups_; };
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
    bool evaluateGradient(const Vector& x, const std::vector<int>& cols,
                          Vector& g);
    /**
     * Evaluates gradient in subspace defined by variables
     * in subgroups.
     * \param[in] x is a given point/iterate, a constant double array
     * \param[in] subgroups is a vector of column indices
     * \param[in] v is a given vector, a constant double array
     * \param[out] Hv is the product of the Hessian and "v", a double array
     * (return value) \return indicator of success (true) or failure (false)
     */
    bool evaluateHessianVectorProduct(const Vector&           x,
                                      const std::vector<int>& cols,
                                      const Vector& v, Vector& Hv);
    /**
     * @brief Evaluator proximal operator at u
     *
     * \param u
     * \param quantities
     * \return true
     * \return false
     */
    bool evaluateProximalGradient(const Vector& x, const Vector& gradfx,
                                  float stepsize, Vector& proxgrad);
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
    float              penalty_;
    std::vector<float> weights_;
    /** pointer to a vector of vector specifing the groups of variables
     *  For example, consider the following group structures:
     *  group1: [0,2,4]; group2: [1,3,5,6]
     *  the groups_ is a pointer to the vector {{0,2,4}, {1,3,5,6}}
     */
    std::vector<std::vector<int>> groups_;
    std::vector<float>            per_group_norm_;
    std::vector<int>              per_group_size_;
    std::map<int, int>            idx_to_group_;
    std::map<int, float>          per_group_gradient_step_norm_;
    //@}

};  // end GroupL1

#endif /* __GROUPL1_HPP__ */
