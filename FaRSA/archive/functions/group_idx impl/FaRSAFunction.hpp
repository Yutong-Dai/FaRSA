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
 * Forward declarations
 */
class Quantities;
class Reporter;

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
    int              numberOfVariables() { return number_of_variables_; };
    inline int const numberOfGroups() const { return number_of_groups_; };
    /**
     * @brief Get a 2-d vector specifies the group structure.
     *
     * \return std::vector<std::vector<int>>&
     */
    inline std::shared_ptr<std::vector<std::vector<int>>> const groups() const
    {
        return groups_;
    };
    std::string const name() const { return name_; };
    //@}

    /** @name Evaluate methods */
    //@{
    /**
     * Evaluates objective
     * \param[in] x is a given point/iterate, a constant double array
     * \param[out] f is the objective value at "x", a double (return value)
     */
    virtual bool evaluateObjective(const Vector& x, double& f) = 0;
    // /**
    //  * Evaluates gradient
    //  * \param[in] x is a given point/iterate, a constant double array
    //  * \param[out] g is the gradient value at "x", a double array (return
    //  value)
    //  */
    // virtual bool evaluateGradient(const Vector& x, Vector& g) = 0;
    /**
     * Evaluates gradient in subspace defined by variables
     * in cols.
     * \param[in] x is a given point/iterate, a constant double array
     * \param[in] group_idx is a vector of column indices; if empty, then
     * calculate full gradient.
     * \param[out] g is the gradient value at "x", a double array (return value)
     */
    virtual bool evaluateGradient(const Vector&           x,
                                  const std::vector<int>& group_idx,
                                  Vector&                 g) = 0;
    /**
     * Evaluates Hessian-vector product in subspace defined by variables
     * in cols.
     * \param[in] x is a given point/iterate, a constant double array
     * \param[in] group_idx is a vector of selected group indices
     * \param[in] v is a given vector, a constant double array
     * \param[out] Hv is the Hessian value at "x" times "v", a double array
     * (return value)
     */
    virtual bool evaluateHessianVectorProduct(const Vector&           x,
                                              const std::vector<int>& group_idx,
                                              const Vector& v, Vector& Hv) = 0;
    //@}

   protected:
    /** @name Protected members */
    //@{
    int         number_of_variables_;
    std::string name_;
    int         number_of_groups_;
    /** a vector of vector specifing the groups of variables
     *  For example, consider the following group structures:
     *  Non-overlapping:
     *      group1: [0,1]; group2: [2,3,4,5]; group3:[6]
     *      the groups_ is set as [[0,1], [2,5], [6,6]].
     *      For i-th subvector,
     *      the first element is the statring position of the i-th group
     *      the second element is the ending position of the i-th group
     *  Overlapping:
     *      group1: [0,1,2]; group2: [2,3,4]; group3:[1,3]
     *      the groups_ is set as [[0,1,2], [2,3,4], [1,3]].
     *  The smooth function's groups_ is determined by the nonmooth function's
     * group_
     */
    std::shared_ptr<std::vector<std::vector<int>>> groups_;
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
