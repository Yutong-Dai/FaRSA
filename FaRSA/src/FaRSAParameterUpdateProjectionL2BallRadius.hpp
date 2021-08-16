// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __FARSAPARAMETERUPDATEPROJECTIONL2BALLRADIUS_HPP__
#define __FARSAPARAMETERUPDATEPROJECTIONL2BALLRADIUS_HPP__

#include <memory>
#include <string>

#include "FaRSAEnumerations.hpp"
#include "FaRSAOptions.hpp"
#include "FaRSAParameterUpdate.hpp"
#include "FaRSAQuantities.hpp"
#include "FaRSAReporter.hpp"
#include "FaRSAStrategies.hpp"
#include "FaRSAStrategy.hpp"

namespace FaRSA
{
/**
 * Forward declarations
 */
class Options;
class Quantities;
class Reporter;
class Strategies;
class Strategy;

/**
 * LineSearch class
 */
class ParameterUpdateProjectionL2BallRadius : public ParameterUpdate
{
   public:
    /** @name Constructors */
    //@{
    /**
     * Constructor
     */
    ParameterUpdateProjectionL2BallRadius() { name_ = "ProjectionL2BallRadius"; };
    //@}

    /** @name Destructor */
    //@{
    /**
     * Destructor
     */
    ~ParameterUpdateProjectionL2BallRadius(){};

    /** @name Options handling methods */
    //@{
    /**
     * Add options
     * \param[in,out] options is pointer to Options object from FaRSA
     * \param[in] reporter is pointer to Reporter object from FaRSA
     */
    void addOptions(Options* options, const Reporter* reporter);
    /**
     * Set options
     * \param[in] options is pointer to Options object from FaRSA
     * \param[in] reporter is pointer to Reporter object from FaRSA
     */
    void getOptions(const Options* options, const Reporter* reporter);
    //@}

    /** @name Initialization method */
    //@{
    /**
     * Initialize strategy
     * \param[in] options is pointer to Options object from FaRSA
     * \param[in,out] quantities is pointer to Quantities object from FaRSA
     * \param[in] reporter is pointer to Reporter object from FaRSA
     */
    void initialize(const Options* options, Quantities* quantities, const Reporter* reporter);
    //@}

    /** @name Get method */
    //@{
    // /**
    //  * Get status
    //  * \return current status of line search strategy
    //  */
    // inline PU_Status status() { return status_; };
    //@}

    /** @name Set method */
    //@{
    // /**
    //  * Set status
    //  * \param[in] status is new status to be set
    //  */
    // inline void setStatus(PU_Status status) { status_ = status; };
    //@}

    /** @name Update method */
    //@{
    /**
     * Update the parameter
     * \param[in] options is pointer to Options object from FaRSA
     * \param[in,out] quantities is pointer to Quantities object from FaRSA
     * \param[in] reporter is pointer to Reporter object from FaRSA
     * \param[in,out] strategies is pointer to Strategies object from FaRSA
     */
    void update(const Options* options, Quantities* quantities, const Reporter* reporter, Strategies* strategies);
    //@}

   private:
    /** @name Default compiler generated methods
     * (Hidden to avoid implicit creation/calling.)
     */
    //@{
    /**
     * Copy constructor
     */
    ParameterUpdateProjectionL2BallRadius(const ParameterUpdateProjectionL2BallRadius&);
    /**
     * Overloaded equals operator
     */
    void operator=(const ParameterUpdateProjectionL2BallRadius&);
    //@}

    /** @name Private members */
    //@{
    double kappa1_max_;
    double kappa1_min_;
    double kappa2_max_;
    double kappa2_min_;
    double kappa_increase_factor_;
    double kappa_decrease_factor_;
    //@}

};  // end ParameterUpdateProjectionL2BallRadius

}  // namespace FaRSA

#endif /* __FARSAPARAMETERUPDATEPROJECTIONL2BALLRADIUS_HPP__ */
