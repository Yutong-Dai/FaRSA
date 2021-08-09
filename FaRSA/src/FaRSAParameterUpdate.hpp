// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __FARSAPARAMETERUPDATE_HPP__
#define __FARSAPARAMETERUPDATE_HPP__

#include <memory>
#include <string>

#include "FaRSAEnumerations.hpp"
#include "FaRSAOptions.hpp"
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
class ParameterUpdate : public Strategy
{
   public:
    /** @name Constructors */
    //@{
    /**
     * Constructor
     */
    ParameterUpdate(){};
    //@}

    /** @name Destructor */
    //@{
    /**
     * Destructor
     */
    virtual ~ParameterUpdate(){};

    /** @name Options handling methods */
    //@{
    /**
     * Add options
     * \param[in,out] options is pointer to Options object from FaRSA
     * \param[in] reporter is pointer to Reporter object from FaRSA
     */
    virtual void addOptions(Options* options, const Reporter* reporter) = 0;
    /**
     * Set options
     * \param[in] options is pointer to Options object from FaRSA
     * \param[in] reporter is pointer to Reporter object from FaRSA
     */
    virtual void getOptions(const Options* options, const Reporter* reporter) = 0;
    //@}

    /** @name Initialization method */
    //@{
    /**
     * Initialize strategy
     * \param[in] options is pointer to Options object from FaRSA
     * \param[in,out] quantities is pointer to Quantities object from FaRSA
     * \param[in] reporter is pointer to Reporter object from FaRSA
     */
    virtual void initialize(const Options* options, Quantities* quantities,
                            const Reporter* reporter) = 0;
    //@}

    /** @name Get method */
    //@{
    /**
     * Get name of strategy
     * \return string with name of strategy
     */
    virtual std::string name() = 0;
    /**
     * Get status
     * \return current status of line search strategy
     */
    inline PU_Status status() { return status_; };

    inline std::string method() { return method_; };
    //@}

    /** @name Set method */
    //@{
    /**
     * Set status
     * \param[in] status is new status to be set
     */
    inline void setStatus(PU_Status status) { status_ = status; };
    inline void setMethod(std::string method) { method_ = method; };
    //@}

    /** @name Update method */
    //@{
    /**
     * Update parameters
     * \param[in] options is pointer to Options object from FaRSA
     * \param[in,out] quantities is pointer to Quantities object from FaRSA
     * \param[in] reporter is pointer to Reporter object from FaRSA
     * \param[in,out] strategies is pointer to Strategies object from FaRSA
     */
    virtual void update(const Options* options, Quantities* quantities, const Reporter* reporter,
                        Strategies* strategies) = 0;
    //@}
   protected:
    /** @name Protected members */
    //@{
    std::string method_;
    //@}
   private:
    /** @name Default compiler generated methods
     * (Hidden to avoid implicit creation/calling.)
     */
    //@{
    /**
     * Copy constructor
     */
    ParameterUpdate(const ParameterUpdate&);
    /**
     * Overloaded equals operator
     */
    void operator=(const ParameterUpdate&);
    //@}

    /** @name Private members */
    //@{
    PU_Status status_;
    //@}

};  // end ParameterUpdate

}  // namespace FaRSA

#endif /* __FARSAPARAMETERUPDATE_HPP__ */
