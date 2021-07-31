// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __FARSAFUNCTION_HPP__
#define __FARSAFUNCTION_HPP__

#include <cmath>
#include <string>
#include <vector>

#include "FaRSAVector.hpp"

namespace FaRSA
{
/**
 * Problem class
 */
class Function
{
   public:
    /** @name Constructors */
    //@{
    /**
     * Constructor
     */
    Function(){};
    //@}

    /** @name Destructor */
    //@{
    /**
     * Destructor
     */
    virtual ~Function(){};

    /** @name Get methods */
    //@{
    int         getNumberOfVariables() { return number_of_variables_; };
    std::string getName() { return name_; };
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
     * \param[out] g is the gradient value at "x", a double array (return value)
     */
    virtual bool evaluateGradient(const Vector& x, Vector& g) = 0;
    /**
     * Evaluates gradient in subspace defined by variables
     * in cols.
     * \param[in] x is a given point/iterate, a constant double array
     * \param[in] cols is a vector of group indices
     * \param[out] g is the gradient value at "x", a double array (return value)
     */
    virtual bool evaluateGradient(const Vector& x, const std::vector<int>& cols,
                                  Vector& g) = 0;
    /**
     * Evaluates Hessian-vector product in subspace defined by variables
     * in cols.
     * \param[in] x is a given point/iterate, a constant double array
     * \param[in] cols is a vector of selected group indices
     * \param[in] v is a given vector, a constant double array
     * \param[out] Hv is the Hessian value at "x" times "v", a double array
     * (return value)
     */
    virtual bool evaluateHessianVectorProduct(const Vector&           x,
                                              const std::vector<int>& cols,
                                              const Vector& v, Vector& Hv) = 0;
    //@}

   protected:
    /** @name Protected members */
    //@{
    int         number_of_variables_;
    std::string name_;
    //@}

   private:
    /** @name Default compiler generated methods
     * (Hidden to avoid implicit creation/calling.)
     */
    //@{
    /**
     * Copy constructor
     */
    Function(const Function&);
    /**
     * Overloaded equals operator
     */
    void operator=(const Function&);
    //@}

    /**
     * Private members
     */
    //@{

    //@}

};  // end Function

}  // namespace FaRSA

#endif /* __FARSAFUNCTION_HPP__ */
