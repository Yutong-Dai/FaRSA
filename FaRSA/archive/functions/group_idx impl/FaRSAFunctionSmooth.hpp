// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __FARSAFUNCTIONSMOOTH_HPP__
#define __FARSAFUNCTIONSMOOTH_HPP__

#include "FaRSAFunction.hpp"
#include "FaRSAMatrix.hpp"
#include "FaRSAVector.hpp"

namespace FaRSA
{
/**
 * Problem class
 */
class FunctionSmooth : public Function
{
   public:
    /** @name Constructors */
    //@{
    /**
     * Constructor
     */
    FunctionSmooth(){};
    //@}

    /** @name Destructor */
    //@{
    /**
     * Destructor
     */
    virtual ~FunctionSmooth(){};

    /** @name Get methods */
    //@{
    //@}

    /** @name Evaluate methods */
    //@{
    //@}

   protected:
    /** @name Protected members */
    //@{
    Matrix      data_matrix_;
    Vector      data_label_;
    std::string dataset_name_;
    int         number_of_data_points_;
    //@}

   private:
    /** @name Default compiler generated methods
     * (Hidden to avoid implicit creation/calling.)
     */
    //@{
    /**
     * Copy constructor
     */
    FunctionSmooth(const FunctionSmooth&);
    /**
     * Overloaded equals operator
     */
    void operator=(const FunctionSmooth&);
    //@}
    /**
     * Private members
     */
    //@{
    //@}

};  // end FunctionSmooth

}  // namespace FaRSA

#endif /* __FARSAFUNCTIONSMOOTH_HPP__ */
