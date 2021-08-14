// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __FARSASPACEPARTITION_HPP__
#define __FARSASPACEPARTITION_HPP__

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
 * SpacePartition class
 */
class SpacePartition : public Strategy
{
   public:
    /** @name Constructors */
    //@{
    /**
     * Constructor
     */
    SpacePartition(){};
    //@}

    /** @name Destructor */
    //@{
    /**
     * Destructor
     */
    virtual ~SpacePartition(){};

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
    virtual void initialize(const Options* options, Quantities* quantities, const Reporter* reporter) = 0;
    //@}

    /** @name Get methods */
    //@{
    /**
     * Get iteration header string
     * \return string of header values
     */
    virtual std::string iterationHeader() = 0;
    /**
     * Get iteration null values string
     * \return string of null values
     */
    virtual std::string iterationNullValues() = 0;
    /**
     * Get name of strategy
     * \return string with name of strategy
     */
    virtual std::string name() = 0;
    /**
     * Get status
     * \return current status of direction computation
     */
    inline SP_Status status() { return status_; };
    //@}

    /** @name Set method */
    //@{
    /**
     * Set status
     * \param[in] status is new status to be set
     */
    inline void setStatus(SP_Status status) { status_ = status; };
    //@}

    /** @name Direction computation method */
    //@{
    /**
     * Run direction computation
     * \param[in] options is pointer to Options object from FaRSA
     * \param[in,out] quantities is pointer to Quantities object from FaRSA
     * \param[in] reporter is pointer to Reporter object from FaRSA
     * \param[in,out] strategies is pointer to Strategies object from FaRSA
     */
    virtual void partitionSpace(const Options* options, Quantities* quantities, const Reporter* reporter,
                                Strategies* strategies) = 0;
    //@}

   private:
    /** @name Default compiler generated methods
     * (Hidden to avoid implicit creation/calling.)
     */
    //@{
    /**
     * Copy constructor
     */
    SpacePartition(const SpacePartition&);
    /**
     * Overloaded equals operator
     */
    void operator=(const SpacePartition&);
    //@}

    /** @name Private members */
    //@{
    SP_Status status_; /**< Termination status */
                       //@}

};  // end SpacePartition

}  // namespace FaRSA

#endif /* __FARSASPACEPARTITION_HPP__ */
