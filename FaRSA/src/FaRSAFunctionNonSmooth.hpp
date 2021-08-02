// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __FARSAFUNCTIONNONSMOOTH_HPP__
#define __FARSAFUNCTIONNONSMOOTH_HPP__

#include "FaRSAFunction.hpp"
#include "FaRSAMatrix.hpp"
#include "FaRSAVector.hpp"

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
    FunctionNonSmooth(){};
    //@}

    /** @name Destructor */
    //@{
    /**
     * Destructor
     */
    virtual ~FunctionNonSmooth(){};

    /** @name Get methods */
    //@{
    //@}

    /** @name Evaluate methods */
    //@{
    // /**
    //  * Evaluates objective
    //  * \param[in] x is a given point/iterate, a constant double array
    //  * \param[out] f is the objective value at "x", a double (return value)
    //  */
    // virtual bool evaluateObjective(const Vector& x, double& f) = 0;
    /**
     * Evaluates gradient. This method is not intended for the nonsmooth
     function, so will always return false.
     * \param[in] x is a given point/iterate, a constant double array
     * \param[out] g is the gradient value at "x", a double array (return
     value)
     */
    bool evaluateGradient(const Vector& x, Vector& g) { return false; };
    /**
     * Evaluates gradient
     * \param[in] x is a given point/iterate, a constant double array
     * \param[in] cols is a vector of group indices
     * \param[out] g is the gradient value at "x", a double array (return
     * value)
     */
    virtual bool evaluateGradient(const Vector& x, const std::vector<int>& cols,
                                  Vector& g) = 0;
    // /**
    //  * Evaluates Hessian-vector product
    //  * \param[in] x is a given point/iterate, a constant double array
    //  * \param[in] cols is a vector of group indices
    //  * \param[in] v is a given vector, a constant double array
    //  * \param[out] Hv is the Hessian value at "x" times "v", a double
    //  array
    //  * (return value)
    //  */
    // virtual bool evaluateHessianVectorProduct(const Vector&          x,
    //                                           const std::vector<int>
    //                                           cols, const Vector& v,
    //                                           Vector& Hv) = 0;

    virtual bool evaluateProximalGradient(const Vector& x, const Vector& gradfx,
                                          double  stepsize,
                                          Vector& proxgrad) = 0;
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

};  // end FunctionNonSmooth

}  // namespace FaRSA

#endif /* __FARSAFUNCTIONNONSMOOTH_HPP__ */
