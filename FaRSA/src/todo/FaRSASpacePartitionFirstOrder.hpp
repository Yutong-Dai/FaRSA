// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __FARSASPACEPARTITIONFIRSTORDER_HPP__
#define __FARSASPACEPARTITIONFIRSTORDER_HPP__

#include "FaRSASpacePartition.hpp"

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
class SpacePartitionFirstOrder : public SpacePartition
{
   public:
    /** @name Constructors */
    //@{
    /**
     * Constructor
     */
    SpacePartitionFirstOrder(){};
    //@}

    /** @name Destructor */
    //@{
    /**
     * Destructor
     */
    ~SpacePartitionFirstOrder(){};

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
    void initialize(const Options* options, Quantities* quantities,
                    const Reporter* reporter);
    //@}

    /** @name Get methods */
    //@{
    /**
     * Get name of strategy
     * \return string with name of strategy
     */
    std::string name() { return "First Order Partition(Full Space)"; };
    //@}

    /** @name Space Partition method */
    //@{
    /**
     * Run direction computation
     * \param[in] options is pointer to Options object from FaRSA
     * \param[in,out] quantities is pointer to Quantities object from FaRSA
     * \param[in] reporter is pointer to Reporter object from FaRSA
     * \param[in,out] strategies is pointer to Strategies object from FaRSA
     */
    void partitionSpace(const Options* options, Quantities* quantities,
                        const Reporter* reporter, Strategies* strategies);
    //@}

   private:
    /** @name Default compiler generated methods
     * (Hidden to avoid implicit creation/calling.)
     */
    //@{
    /**
     * Copy constructor
     */
    SpacePartitionFirstOrder(const SpacePartitionFirstOrder&);
    /**
     * Overloaded equals operator
     */
    void operator=(const SpacePartitionFirstOrder&);
    //@}

};  // end SpacePartition

}  // namespace FaRSA

#endif /* __FARSASPACEPARTITIONFIRSTORDER_HPP__ */
