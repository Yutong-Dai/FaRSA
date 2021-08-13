// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __FARSAPARAMETERUPDATEPGSTEPSIZE_HPP__
#define __FARSAPARAMETERUPDATEPGSTEPSIZE_HPP__

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
class ParameterUpdatePGStepsize : public ParameterUpdate
{
   public:
    /** @name Constructors */
    //@{
    /**
     * Constructor
     */
    ParameterUpdatePGStepsize() { name_ = "PGStepsize"; };
    //@}

    /** @name Destructor */
    //@{
    /**
     * Destructor
     */
    ~ParameterUpdatePGStepsize(){};

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
    ParameterUpdatePGStepsize(const ParameterUpdatePGStepsize&);
    /**
     * Overloaded equals operator
     */
    void operator=(const ParameterUpdatePGStepsize&);
    //@}

    /** @name Private members */
    //@{
    // PU_Status status_;
    double decrease_factor_;
    double increase_factor_;
    double upper_bound_;
    //@}

};  // end ParameterUpdatePGStepsize

}  // namespace FaRSA

#endif /* __FARSAPARAMETERUPDATEPGSTEPSIZE_HPP__ */
