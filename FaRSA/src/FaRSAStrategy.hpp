// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __FARSASTRATEGY_HPP__
#define __FARSASTRATEGY_HPP__

#include <string>

#include "FaRSAOptions.hpp"
#include "FaRSAQuantities.hpp"
#include "FaRSAReporter.hpp"

namespace FaRSA
{
/**
 * Forward declarations
 */
class Options;
class Quantities;
class Reporter;

/**
 * Strategy class
 */
class Strategy
{
   public:
    /** @name Constructors */
    //@{
    /**
     * Constructor
     */
    Strategy(){};
    //@}

    /** @name Destructor */
    //@{
    /**
     * Destructor
     */
    virtual ~Strategy(){};

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

    /** @name Get method */
    //@{
    /**
     * Get name of strategy
     * \return string with name of strategy
     */
    virtual std::string name() = 0;
    //@}

   private:
    /** @name Default compiler generated methods
     * (Hidden to avoid implicit creation/calling.)
     */
    //@{
    /**
     * Copy constructor
     */
    Strategy(const Strategy&);
    /**
     * Overloaded equals operator
     */
    void operator=(const Strategy&);
    //@}

};  // end Strategy

}  // namespace FaRSA

#endif /* __FARSASTRATEGY_HPP__ */
