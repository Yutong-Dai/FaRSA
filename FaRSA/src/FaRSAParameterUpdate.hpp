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
 * ParameterUpdate class
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
    virtual void initialize(const Options* options, Quantities* quantities, const Reporter* reporter) = 0;
    //@}

    /** @name Get method */
    //@{
    /**
     * Get name of strategy
     * \return string with name of strategy
     */
    inline std::string name() { return name_; };
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
    std::string name_;
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

/**
 * ParameterUpdates class: host a list of ParameterUpdate classs
 */
class ParameterUpdates
{
   public:
    /** @name Constructors */
    //@{
    /**
     * Constructor
     */
    ParameterUpdates(){};
    //@}

    /** @name Destructor */
    //@{
    /**
     * Destructor
     */
    ~ParameterUpdates(){};
    //@}

    /** @name Add methods */
    //@{
    /**
     * Add bool option
     * \param[in] reporter is pointer to Reporter object from FaRSA
     * \param[in] parameter_update is a pointer to ParameterUpdate class
     * \param[in] description is description of option
     */
    bool add(std::shared_ptr<ParameterUpdate> parameter_update)
    {
        int old_len = parameter_update_list_.size();
        parameter_update_list_.push_back(parameter_update);
        if (parameter_update_list_.size() == old_len + 1)
        {
            details_ += "|-- " + parameter_update->name() + ": " + parameter_update->method() + "\n";
            return true;
        }
        else
        {
            return false;
        }
    };
    //@}

    /** @name Get method */
    //@{
    /**
     * Get name of a set of parameter update strategies
     * \return string with name of parameter update strategies
     */
    inline const std::string details() const { return details_; };
    inline const PU_Status   statusOverall() const { return status_overall_; };
    inline const int         size() { return parameter_update_list_.size(); };
    //@}

    /** @name Initialization method */
    //@{
    /**
     * Initialize strategy
     * \param[in] options is pointer to Options object from FaRSA
     * \param[in,out] quantities is pointer to Quantities object from FaRSA
     * \param[in] reporter is pointer to Reporter object from FaRSA
     */
    void initialize(const Options* options, Quantities* quantities, const Reporter* reporter)
    {
        for (auto parameter_update : parameter_update_list_)
        {
            parameter_update->initialize(options, quantities, reporter);
        }
    }
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
    inline void update(const Options* options, Quantities* quantities, const Reporter* reporter, Strategies* strategies)
    {
        for (auto parameter_update : parameter_update_list_)
        {
            parameter_update->update(options, quantities, reporter, strategies);
            if (parameter_update->status() != PU_SUCCESS)
            {
                status_overall_ = PU_UPDATE_FAILURE;
                reporter->printf(R_SOLVER, R_BASIC, "Parameter update strategy %s failed.\n",
                                 parameter_update->name().c_str());
            }
        }
        status_overall_ = PU_SUCCESS;
    }
    //@}

    //    protected:
    //     /** @name Protected members */
    //     //@{

    //     //@}
   private:
    /** @name Default compiler generated methods
     * (Hidden to avoid implicit creation/calling.)
     */
    //@{
    /**
     * Copy constructor
     */
    ParameterUpdates(const ParameterUpdates&);
    /**
     * Overloaded equals operator
     */
    void operator=(const ParameterUpdates&);
    //@}

    /** @name Private members */
    //@{
    std::vector<std::shared_ptr<ParameterUpdate>>
                parameter_update_list_; /**< Vector of (pointers to) ParameterUpdate */
    std::string details_;
    PU_Status   status_overall_;
    //@}

};  // end ParameterUpdates

}  // namespace FaRSA

#endif /* __FARSAPARAMETERUPDATE_HPP__ */
