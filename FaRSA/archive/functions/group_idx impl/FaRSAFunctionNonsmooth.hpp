// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __FARSAFunctionNonsmooth_HPP__
#define __FARSAFunctionNonsmooth_HPP__

#include "FaRSAFunction.hpp"
#include "FaRSAMatrix.hpp"
#include "FaRSAQuantities.hpp"
#include "FaRSAVector.hpp"

namespace FaRSA
{
class FunctionNonsmooth : public Function
{
   public:
    /** @name Constructors */
    //@{
    /**
     * Constructor
     */
    FunctionNonsmooth(){};
    //@}

    /** @name Destructor */
    //@{
    /**
     * Destructor
     */
    virtual ~FunctionNonsmooth(){};

    // /** @name Get methods */
    // //@{

    // //@}

    /** @name Evaluate methods */
    //@{
    //@}

    /** @name Compute methods */
    //@{
    virtual bool computeProximalGradientUpdate(Quantities& quantities) = 0;
    //@}

    //    protected:
    //     /** @name Protected members */
    //     //@{

    //     //@}

   private:
    /** @name Default compiler generated methods
     * (Hidden to avoid implicit creation/calling.)
     */
    //@{
    /**
     * Copy constructor
     */
    FunctionNonsmooth(const FunctionNonsmooth&);
    /**
     * Overloaded equals operator
     */
    void operator=(const FunctionNonsmooth&);
    //@}

};  // end FunctionNonsmooth

}  // namespace FaRSA

#endif /* __FARSAFunctionNonsmooth_HPP__ */
