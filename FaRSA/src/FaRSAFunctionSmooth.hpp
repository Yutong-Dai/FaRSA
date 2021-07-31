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
     * Evaluates gradient. This method is not intended for the smooth function,
     * so will always return false.
     * \param[in] x is a given point/iterate, a constant double array
     * \param[in] cols is a vector of group indices
     * \param[out] g is the gradient value at "x", a double array (return value)
     */
    bool evaluateGradient(const Vector& x, const std::vector<int>& cols,
                          Vector& g)
    {
        return false;
    };
    /**
     * Evaluates Hessian-vector product in subspace defined by variables
     * in cols.
     * \param[in] x is a given point/iterate, a constant double array
     * \param[in] cols is a vector of group indices
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
