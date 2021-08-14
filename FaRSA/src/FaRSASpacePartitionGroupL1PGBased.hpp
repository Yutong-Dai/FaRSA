// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __FARSASPACEPARTITIONGROUPL1PGBASED_HPP__
#define __FARSASPACEPARTITIONGROUPL1PGBASED_HPP__

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
 * SpacePartitionGroupL1PGBased class
 */
class SpacePartitionGroupL1PGBased : public SpacePartition
{
   public:
    /** @name Constructors */
    //@{
    /**
     * Constructor
     */
    SpacePartitionGroupL1PGBased()
        : p_(FARSA_DOUBLE_INFINITY),
          kappa1_(FARSA_DOUBLE_INFINITY),
          kappa2_(FARSA_DOUBLE_INFINITY),
          gamma_(FARSA_DOUBLE_INFINITY),
          verbose_(true){};
    //@}

    /** @name Destructor */
    //@{
    /**
     * Destructor
     */
    ~SpacePartitionGroupL1PGBased(){};

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
    std::string iterationHeader()
    {
        if (verbose_)
        {
            return "ProxStepsize #1stGrps #2ndGrps 1stOptim. 2ndOptim.";
        }
        else
        {
            return "ProxStepsize 1stOptim. 2ndOptim.";
        }
    };
    /**
     * Get iteration null values string
     * \return string of null values
     */
    std::string iterationNullValues()
    {
        if (verbose_)
        {
            return " ---------  --------- --------- --------- ---------";
        }
        else
        {
            return " ---------  --------- ---------";
        };
    }
    /**
     * Get name of strategy
     * \return string with name of strategy
     */
    std::string name() { return "GroupL1PGBasedPartition"; };
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
    void partitionSpace(const Options* options, Quantities* quantities, const Reporter* reporter,
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
    SpacePartitionGroupL1PGBased(const SpacePartitionGroupL1PGBased&);
    /**
     * Overloaded equals operator
     */
    void operator=(const SpacePartitionGroupL1PGBased&);
    //@}
    /** @name Private members */
    //@{
    // parameter controls the global iteation complexity result
    double                            p_;
    double                            kappa1_;
    double                            kappa2_;
    double                            gamma_;
    std::shared_ptr<std::vector<int>> groups_full_indices_;
    std::shared_ptr<std::vector<int>> variables_full_indices_;
    bool                              verbose_;

    //@}
};  // end SpacePartitionGroupL1PGBased

}  // namespace FaRSA

#endif /* __FARSASPACEPARTITIONGROUPL1PGBASED_HPP__ */
