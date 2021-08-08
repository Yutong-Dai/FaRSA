// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __FARSADIRECTIONCOMPUTATION_HPP__
#define __FARSADIRECTIONCOMPUTATION_HPP__

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
 * DirectionComputation class
 */
class DirectionComputation : public Strategy
{
   public:
    /** @name Constructors */
    //@{
    /**
     * Constructor
     */
    DirectionComputation(){};
    //@}

    /** @name Destructor */
    //@{
    /**
     * Destructor
     */
    virtual ~DirectionComputation(){};

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
    inline DC_Status status() { return status_; };
    inline bool      performCompuation() { return perform_computation_; };
    //@}

    /** @name Set method */
    //@{
    /**
     * Set status
     * \param[in] status is new status to be set
     */
    inline void setStatus(DC_Status status) { status_ = status; };

    void setPerformComputation(bool perform_computation)
    {
        perform_computation_ = perform_computation;
    }
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
    virtual void computeDirection(const Options* options, Quantities* quantities,
                                  const Reporter* reporter, Strategies* strategies) = 0;
    //@}

    //    protected:
    /** @name Protected members */
    //@{

    //@}

   private:
    /** @name Default compiler generated methods
     * (Hidden to avoid implicit creation/calling.)
     */
    //@{
    /**
     * Copy constructor
     */
    DirectionComputation(const DirectionComputation&);
    /**
     * Overloaded equals operator
     */
    void operator=(const DirectionComputation&);
    //@}

    /** @name Private members */
    //@{
    DC_Status status_; /**< Termination status */
    bool      perform_computation_;
    //@}

};  // end DirectionComputation

}  // namespace FaRSA

#endif /* __FARSADIRECTIONCOMPUTATION_HPP__ */
