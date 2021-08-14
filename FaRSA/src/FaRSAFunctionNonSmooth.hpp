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
    inline const int numberOfGroups() const { return number_of_groups_; };

    inline const double penalty() const { return penalty_; };
    /**
     * @brief Get a 2-d vector specifies the group structure.
     *
     * \return std::vector<std::vector<int>>&
     */
    inline const std::shared_ptr<std::vector<std::vector<int>>> groups() const { return groups_; };
    virtual const std::string                                   groupsFormat() const = 0;
    // //@}

    /** @name Evaluate methods */
    //@{
    //@}
    /** @name scale methods */
    //@{
    virtual bool scalePenalty(double scale) = 0;
    //@}

    /** @name Compute methods */
    //@{
    virtual bool computeProximalGradientUpdate(const Vector& x, const Vector& gradfx, Quantities& quantities,
                                               Vector& proxgrad) = 0;
    //@}

   protected:
    /** @name Protected members */
    //@{
    int number_of_groups_;
    /** a vector of vector specifing the groups of variables
     *  For example, consider the following group structures:
     *  Non-overlapping:
     *      group1: [0,1]; group2: [2,3,4,5]; group3:[6]
     *      the groups_ is set as [[0,1], [2,5], [6,6]].
     *      For i-th subvector,
     *      the first element is the statring position of the i-th group
     *      the second element is the ending position of the i-th group
     *  Overlapping:
     *      group1: [0,1,3]; group2: [2,3,4]; group3:[0,2,4]
     *      the groups_ is set as [[0,1,3], [2,3,4], [0,2,4]].
     * Within each group, the indicies have to be ordered. The code doesn't check this.
     */
    std::shared_ptr<std::vector<std::vector<int>>> groups_;
    double                                         penalty_;
    //@}

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
