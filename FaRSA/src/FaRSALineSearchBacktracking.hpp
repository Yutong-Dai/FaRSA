// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __FARSALINESEARCHBACKTRACKING_HPP__
#define __FARSALINESEARCHBACKTRACKING_HPP__

#include "FaRSALineSearch.hpp"

namespace FaRSA
{
/**
 * LineSearchBacktracking class
 */
class LineSearchBacktracking : public LineSearch
{
   public:
    /** @name Constructors */
    //@{
    /**
     * Constructor
     */
    LineSearchBacktracking(){};
    //@}

    /** @name Destructor */
    //@{
    /**
     * Destructor
     */
    ~LineSearchBacktracking(){};

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
     * Get iteration header values
     * \return string of header values
     */
    std::string iterationHeader()
    {
        if (verbose_)
        {
            if (method_pure_second_order_.compare("projected_armijo_groupl1") == 0)
            {
                return " type #newGrps |Direction| LSStepsize";
            }
            else
            {
                return " |Direction| LSStepsize";
            }
        }
        else
        {
            return " |Direction| LSStepsize";
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
            if (method_pure_second_order_.compare("projected_armijo_groupl1") == 0)
            {
                return " ---- ----- --------- ---------";
            }
            else
            {
                return "--------- ---------";
            }
        }
        else
        {
            return "--------- ---------";
        }
    };
    /**
     * Get name of strategy
     * \return string with name of strategy
     */
    std::string name() { return "Backtracking"; };

    inline int const numberOfBacktrack() const { return number_of_backtrack_; };
    //@}

    /** @name set methods */
    //@{
    inline void reset() { number_of_backtrack_ = 0; };
    //@}

    /** @name Line search method */
    //@{
    /**
     * Run line search
     * \param[in] options is pointer to Options object from FaRSA
     * \param[in,out] quantities is pointer to Quantities object from FaRSA
     * \param[in] reporter is pointer to Reporter object from FaRSA
     * \param[in,out] strategies is pointer to Strategies object from FaRSA
     */
    void runLineSearch(const Options* options, Quantities* quantities, const Reporter* reporter,
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
    LineSearchBacktracking(const LineSearchBacktracking&);
    /**
     * Overloaded equals operator
     */
    void operator=(const LineSearchBacktracking&);
    //@}

    /** @name Private members */
    //@{
    bool        fail_on_small_stepsize_;
    bool        verbose_;
    double      stepsize_initial_;
    double      stepsize_minimum_;
    double      stepsize_sufficient_decrease_threshold_;
    double      stepsize_sufficient_decrease_fudge_factor_;
    double      stepsize_decrease_factor_;
    std::string direction_search_type_;
    std::string directional_derivative_first_order_type_;
    std::string directional_derivative_second_order_type_;
    std::string method_pure_first_order_;
    std::string method_pure_second_order_;
    std::string method_hybrid_;
    // int    max_backtrack_;
    int number_of_backtrack_;
    //@}

};  // end LineSearchBacktracking

}  // namespace FaRSA

#endif /* __FARSALINESEARCHBACKTRACKING_HPP__ */
