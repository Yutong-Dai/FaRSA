// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __FARSADIRECTIONCOMPUTATIONPROXIMALGRADIENT_HPP__
#define __FARSADIRECTIONCOMPUTATIONPROXIMALGRADIENT_HPP__

#include "FaRSADirectionComputation.hpp"

namespace FaRSA
{
/**
 * DirectionComputationProximalGradient class
 */
class DirectionComputationProximalGradient : public DirectionComputation
{
   public:
    /** @name Constructors */
    //@{
    /**
     * Constructor
     */
    DirectionComputationProximalGradient(){};
    //@}

    /** @name Destructor */
    //@{
    /**
     * Destructor
     */
    ~DirectionComputationProximalGradient(){};

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

    /** @name Get methods */
    //@{
    /**
     * Get iteration header string
     * \return string of header values
     */
    std::string iterationHeader();
    /**
     * Get iteration null values string
     * \return string of null values
     */
    std::string iterationNullValues();
    /**
     * Get name of strategy
     * \return string with name of strategy
     */
    std::string name() { return "Proximal Gradient"; };
    //@}

    /** @name Set methods */
    //@{

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
    void computeDirection(const Options* options, Quantities* quantities, const Reporter* reporter,
                          Strategies* strategies);
    //@}

   private:
    /** @name Default compiler generated methods
     * (Hidden to avoid implicit creation/calling.)
     */
    //@{
    /**
     * Copy constructor
     */
    DirectionComputationProximalGradient(const DirectionComputationProximalGradient&);
    /**
     * Overloaded equals operator
     */
    void operator=(const DirectionComputationProximalGradient&);
    //@}

    /** @name Private members */
    //@{
    //@}

};  // end DirectionComputationProximalGradient

}  // namespace FaRSA

#endif /* __FARSADIRECTIONCOMPUTATIONPROXIMALGRADIENT_HPP__ */
