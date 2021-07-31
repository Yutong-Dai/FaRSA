// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __FARSAFUNCTIONNONSMOOTH_HPP__
#define __FARSAFUNCTIONNONSMOOTH_HPP__

#include <cmath>

#include "FaRSAFunction.hpp"
#include "FaRSAPoint.hpp"

namespace FaRSA
{
/**
 * Problem class
 */
class FunctionNonSmooth : public Function
{
   public:
    /** @name Constructors */
    //@{
    /**
     * Constructor
     */
    FunctionNonSmooth(std::vector<std::vector<int>>* groups, float penalty)
    {
        groups_ = groups;
        penalty_ = penalty;
        weights_ = &std::vector<float>(groups_->size());
        for (int i = 0; i < groups_->size(); i++)
        {
            // weight_i = sqrt{|group_i|} * penalty
            float number_of_elements = (float)((*groups_)[i].size());
            (*weights_)[i] = std::sqrt(number_of_elements) * penalty_;
        }
    };
    //@}

    /** @name Destructor */
    //@{
    /**
     * Destructor
     */
    virtual ~FunctionNonSmooth(){};

    /** @name Get methods */
    //@{
    float                          getPenalty() { return penalty_; };
    std::vector<float>*            getWeights() { return weights_; };
    std::vector<std::vector<int>>* getGroups() { return groups_; };
    //@}

    /** @name Evaluate methods */
    //@{
    /**
     * Evaluates objective
     * \param[in] x is a given point/iterate, a constant double array
     * \param[out] f is the objective value at "x", a double (return value)
     */
    virtual bool evaluateObjective(const Vector& x, double& f) = 0;
    /**
     * Evaluates gradient
     * \param[in] x is a given point/iterate, a constant double array
     * \param[in] cols is a vector of group indices
     * \param[out] g is the gradient value at "x", a double array (return value)
     */
    virtual bool evaluateGradient(const Vector& x, const std::vector<int> cols,
                                  Vector& g) = 0;
    /**
     * Evaluates Hessian-vector product
     * \param[in] x is a given point/iterate, a constant double array
     * \param[in] cols is a vector of group indices
     * \param[in] v is a given vector, a constant double array
     * \param[out] Hv is the Hessian value at "x" times "v", a double array
     * (return value)
     */
    virtual bool evaluateHessianVectorProduct(const Vector&          x,
                                              const std::vector<int> cols,
                                              const Vector& v, Vector& Hv) = 0;
    virtual bool evaluateProximalGradient(const std::shared_ptr<Point> point,
                                          float stepsize);
    //@}

    //  protected:
    //   /** @name Protected members */
    //   //@{
    //   int number_of_variables_;               /**< Number of variables */
    //   std::vector<std::vector<int> > groups_; /**< Group data          */
    //                                           //@}

   private:
    /** @name Default compiler generated methods
     * (Hidden to avoid implicit creation/calling.)
     */
    //@{
    /**
     * Copy constructor
     */
    FunctionNonSmooth(const FunctionNonSmooth&);
    /**
     * Overloaded equals operator
     */
    void operator=(const FunctionNonSmooth&);
    //@}
    /**
     * Private members
     */
    //@{
    float               penalty_;
    std::vector<float>* weights_;
    /** pointer to a vector of vector specifing the groups of variables
     *  For example, consider the following group structures:
     *  group1: [0,2,4]; group2: [1,3,5,6]
     *  the groups_ is a pointer to the vector {{0,2,4}, {1,3,5,6}}
     */
    std::vector<std::vector<int>>* groups_;

    //@}

};  // end FunctionNonSmooth

}  // namespace FaRSA

#endif /* __FARSAFUNCTIONNONSMOOTH_HPP__ */
