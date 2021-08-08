// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __FARSAITERATIONQUANTITIES_HPP__
#define __FARSAITERATIONQUANTITIES_HPP__

#include <ctime>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "FaRSAFunctionNonsmooth.hpp"
#include "FaRSAFunctionSmooth.hpp"
#include "FaRSAOptions.hpp"
#include "FaRSAPoint.hpp"
#include "FaRSAReporter.hpp"
#include "FaRSAVector.hpp"

namespace FaRSA
{
/**
 * Forward declarations
 */
class Options;
class Point;
class Reporter;
class Vector;
class FunctionSmooth;
class FunctionNonsmooth;
/**
 * Quantities class
 */
class Quantities
{
   public:
    /** @name Constructors */
    //@{
    /**
     * Declare Quantities
     */
    Quantities();
    //@}

    /** @name Destructor */
    //@{
    /**
     * Delete data
     */
    ~Quantities();
    //@}

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

    /** @name Initialization methods */
    //@{
    /**
     * Initialize quantities
     * \param[in] problem is pointer to Problem object
     */
    bool initialize(const std::shared_ptr<FunctionSmooth>    function_smooth,
                    const std::shared_ptr<FunctionNonsmooth> function_nonsmooth,
                    const Vector&                            initail_point);
    //@}

    /** @name Get methods */
    //@{
    /**
     * Start time
     * \return start time that was set
     */
    inline clock_t const startTime() const { return start_time_; };
    /**
     * \return end time that was set
     */
    inline clock_t const endTime() const { return end_time_; };
    /**
     * Evaluation time
     * \return problem function evaluation time that was set
     */
    inline clock_t const evaluationTime() const { return evaluation_time_; };
    /**
     * Get CPU time limit
     * \return CPU time limit
     */
    inline double const cpuTimeLimit() const { return cpu_time_limit_; };
    /**
     * Get scaling threshold
     * \return scaling threshold
     */
    inline double const scalingThreshold() const { return scaling_threshold_; };
    /**
     * Get stepsize
     * \return current stepsize
     */
    inline double const stepsizeLineSearch() const { return stepsize_ls_; };
    inline double const stepsizeProximalGradient() const { return stepsize_prox_; };
    /**
     * Function evaluation counter
     * \return function evaluations performed so far
     */
    inline int const functionCounter() const { return function_counter_; };
    /**
     * Function evaluation limit
     * \return function evaluation limit
     */
    inline int const functionEvaluationLimit() const { return function_evaluation_limit_; };
    /**
     * Gradient evaluation counter
     * \return gradient evaluations performed so far
     */
    inline int const gradientCounter() const { return gradient_counter_; };

    inline int const hessianVectorCounter() const { return hessian_vector_counter_; };
    /**
     * Gradient evaluation limit
     * \return gradient evaluation limit
     */
    inline int const gradientEvaluationLimit() const { return gradient_evaluation_limit_; };

    inline int const hessianVectorProductEvaluationLimit() const
    {
        return hessian_vector_product_evaluation_limit_;
    };
    /**
     * Iteration counter
     * \return iterations performed so far
     */
    inline int const iterationCounter() const { return iteration_counter_; };
    /**
     * Get problem size
     * \return number of variables
     */
    inline int const numberOfVariables() const { return number_of_variables_; };
    /**
     * Get pointer to current iterate
     * \return pointer to Point representing current iterate
     */
    inline std::shared_ptr<Point> currentIterate() { return current_iterate_; };
    /**
     * Get pointer to trial iterate
     * \return pointer to Point representing trial iterate
     */
    inline std::shared_ptr<Point> trialIterate() { return trial_iterate_; };
    /**
     * Get direction
     * \return pointer to Vector representing search direction
     */
    inline std::shared_ptr<Vector> direction() { return direction_; };

    /**
     * @brief Get stationarity tolerance
     * \return double const stationarity_tolerance_
     */
    inline double const stationarityTolerance() const { return stationarity_tolerance_; };
    /**
     * @brief Get stationarity iteration limit
     *
     * \return double const iteration_limit_
     */
    inline double const iterationLimit() const { return iteration_limit_; };
    /**
     * @brief Get iterate norm tolerance
     *
     * \return double const iterate_norm_tolerance_
     */
    inline double const iterateNormTolerance() const { return iterate_norm_tolerance_; };

    inline double const kappa1Max() const { return kappa1_max_; };
    inline double const kappa1min() const { return kappa1_min_; };
    inline double const kappa2Max() const { return kappa2_max_; };
    inline double const kappa2min() const { return kappa2_min_; };
    inline double const kappaIncreaseFactor() const { return kappa_increase_factor_; };
    inline double const kappaDecreaseFactor() const { return kappa_decrease_factor_; };

    inline std::shared_ptr<std::vector<int>> const groupsFirstOrder() const
    {
        return groups_first_order_;
    };

    inline std::shared_ptr<std::vector<int>> const groupsSecondOrder() const
    {
        return groups_second_order_;
    };

    inline std::shared_ptr<std::vector<int>> const groupsWorking() const
    {
        return groups_working_;
    };
    inline std::shared_ptr<std::vector<int>> const indiciesWorking() const
    {
        return indicies_working_;
    };

    inline int const numberOfGroups() const { return number_of_groups_; };

    inline float const scaleApplied() const { return scale_applied_; };
    inline DC_Type     directionType() { return direction_type_; };
    //@}

    /** @name Set methods */
    //@{
    /**
     * Set current iterate pointer
     * \param[in] iterate is pointer to Point to represent current iterate
     */
    inline void setCurrentIterate(const std::shared_ptr<Point> iterate)
    {
        current_iterate_ = iterate;
    };
    /**
     * Set trial iterate pointer
     * \param[in] trial_iterate is pointer to Point to represent trial iterate
     */
    inline void setTrialIterate(const std::shared_ptr<Point> trial_iterate)
    {
        trial_iterate_ = trial_iterate;
    };
    /**
     * Set trial iterate pointer to current iterate pointer
     */
    inline void setTrialIterateToCurrentIterate() { trial_iterate_ = current_iterate_; };
    /**
     * Set stepsize
     * \param[in] stepsize is new value to represent stepsize
     */
    inline void setStepsizeLineSearch(double stepsize) { stepsize_ls_ = stepsize; };

    inline void setStepsizeProximalGradient(double stepsize) { stepsize_prox_ = stepsize; };

    inline void setGroupsFirstOrder(const std::shared_ptr<std::vector<int>> groups_first_order)
    {
        groups_first_order_ = groups_first_order;
    };

    inline void setGroupsSecondOrder(const std::shared_ptr<std::vector<int>> groups_second_order)
    {
        groups_second_order_ = groups_second_order;
    };

    inline void setGroupsWorking(std::shared_ptr<std::vector<int>> groups_working)
    {
        groups_working_ = groups_working;
    };

    inline void setIndiciesWorking(std::shared_ptr<std::vector<int>> indicies_working)
    {
        indicies_working_ = indicies_working;
    };

    inline void setNumberOfVariables(int number_of_variables)
    {
        number_of_variables_ = number_of_variables;
    };

    inline void setNumberOfGroups(int number_of_groups) { number_of_groups_ = number_of_groups; };

    inline void setScalingThreshold(double scale) { scaling_threshold_ = scale; };
    inline void setScaleApplied(double scale) { scale_applied_ = scale; };
    inline void setDirectionType(DC_Type direction_type) { direction_type_ = direction_type; };
    //@}

    /** @name Increment methods */
    //@{
    /**
     * Increment evaluation time
     * \param[in] evaluation_time is amount to add to total problem function
     * evaluation time
     */
    inline void incrementEvaluationTime(clock_t evaluation_time)
    {
        evaluation_time_ += evaluation_time;
    };
    /**
     * Increment function evaluation counter
     */
    inline void incrementFunctionCounter() { function_counter_++; };
    /**
     * Increment gradient evaluation counter
     */
    inline void incrementGradientCounter() { gradient_counter_++; };
    inline void incrementHessianVectorCounter() { hessian_vector_counter_++; };
    /**
     * Increment iteration counter
     */
    inline void incrementIterationCounter() { iteration_counter_++; };
    //@}

    /** @name Print methods */
    //@{
    /**
     * @brief print values of private members that are loaded from the Option
     * instance.
     *
     * \param[in] reporter is pointer to Reporter object from FaRSA
     */
    void print(const Reporter* reporter);
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
    //@}

    /** @name Print methods */
    //@{
    /**
     * Print header
     * \param[in] reporter is pointer to Reporter object from FaRSA
     */
    void printHeader(const Reporter* reporter);
    /**
     * Print iteration values
     * \param[in] reporter is pointer to Reporter object from FaRSA
     */
    void printIterationValues(const Reporter* reporter);
    /**
     * Print footer
     * \param[in] reporter is pointer to Reporter object from FaRSA
     */
    void printFooter(const Reporter* reporter);
    //@}

    /** @name Finalization method */
    //@{
    /**
     * Finalize
     */
    void finalize();
    //@}

   private:
    /** @name Default compiler generated methods
     * (Hidden to avoid implicit creation/calling.)
     */
    //@{
    /**
     * Copy constructor
     */
    Quantities(const Quantities&);
    /**
     * Overloaded equals operator
     */
    void operator=(const Quantities&);
    //@}

    /** @name Private members */
    //@{
    // generic members
    clock_t                                        end_time_;
    clock_t                                        evaluation_time_;
    clock_t                                        start_time_;
    double                                         stepsize_ls_;
    double                                         stepsize_prox_;
    double                                         scale_applied_;
    int                                            function_counter_;
    int                                            gradient_counter_;
    int                                            hessian_vector_counter_;
    int                                            iteration_counter_;
    int                                            number_of_variables_;
    int                                            number_of_groups_;
    std::shared_ptr<Point>                         current_iterate_;
    std::shared_ptr<Point>                         trial_iterate_;
    std::shared_ptr<Vector>                        direction_;
    std::shared_ptr<std::vector<int>>              groups_first_order_;
    std::shared_ptr<std::vector<int>>              groups_second_order_;
    std::shared_ptr<std::vector<int>>              groups_working_;
    std::shared_ptr<std::vector<int>>              indicies_working_;
    std::shared_ptr<std::vector<std::vector<int>>> groups_;
    DC_Type                                        direction_type_;

    /** @name Private members (set by options) */
    //@{
    double scaling_threshold_;
    // for termination
    double cpu_time_limit_;
    double iterate_norm_tolerance_;
    double stationarity_tolerance_;
    int    function_evaluation_limit_;
    int    gradient_evaluation_limit_;
    int    hessian_vector_product_evaluation_limit_;
    int    iteration_limit_;
    // for space partition using (FaRSAGroup)
    double kappa1_max_;
    double kappa1_min_;
    double kappa2_max_;
    double kappa2_min_;
    double kappa_increase_factor_;
    double kappa_decrease_factor_;
    //@}

};  // end Quantities

}  // namespace FaRSA

#endif /* __FARSAITERATIONQUANTITIES_HPP__ */
