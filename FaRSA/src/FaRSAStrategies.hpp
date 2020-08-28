// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __FARSASTRATEGIES_HPP__
#define __FARSASTRATEGIES_HPP__

#include <memory>
#include <string>

#include "FaRSADirectionComputation.hpp"
#include "FaRSALineSearch.hpp"

namespace FaRSA
{

/**
 * Forward declarations
 */
class DirectionComputation;
class LineSearch;

/**
 * Strategies class
 */
class Strategies
{

public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   */
  Strategies(){};
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  ~Strategies(){};
  //@}

  /** @name Options handling methods */
  //@{
  /**
   * Add options
   * \param[in,out] options is pointer to Options object from FaRSA
   * \param[in] reporter is pointer to Reporter object from FaRSA
   */
  void addOptions(Options* options,
                  const Reporter* reporter);
  /**
   * Get options
   * \param[in] options is pointer to Options object from FaRSA
   * \param[in] reporter is pointer to Reporter object from FaRSA
   */
  void getOptions(const Options* options,
                  const Reporter* reporter);
  //@}

  /** @name Initialization method */
  //@{
  /**
   * Initialize strategies
   * \param[in] options is pointer to Options object from FaRSA
   * \param[in,out] quantities is pointer to Quantities object from FaRSA
   * \param[in] reporter is pointer to Reporter object from FaRSA
   */
  void initialize(const Options* options,
                  Quantities* quantities,
                  const Reporter* reporter);
  //@}

  /** @name Get methods */
  //@{
  /**
   * Get pointer to DirectionComputation
   * \return pointer to DirectionComputation object
   */
  inline std::shared_ptr<DirectionComputation> directionComputation() { return direction_computation_; }
  /**
   * Get pointer to LineSearch
   * \return pointer to LineSearch object
   */
  inline std::shared_ptr<LineSearch> lineSearch() { return line_search_; }
  /**
   * Get iteration header
   * \return iteration header as string
   */
  inline std::string iterationHeader() const { return iteration_header_; }
  //@}

  /** @name Set methods */
  //@{
  /**
   * Set iteration header
   */
  void setIterationHeader();
  //@}

  /** @name Print methods */
  //@{
  /**
   * Print header
   * \param[in] reporter is pointer to Reporter object from FaRSA
   */
  void printHeader(const Reporter* reporter);
  /**
   * Print footer
   * \param[in] reporter is pointer to Reporter object from FaRSA
   */
  void printFooter(const Reporter* reporter);
  //@}

private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Copy constructor
   */
  Strategies(const Strategies&);
  /**
   * Overloaded equals operator
   */
  void operator=(const Strategies&);
  //@}

  /** @name * Private members */
  //@{
  std::shared_ptr<DirectionComputation> direction_computation_;
  std::shared_ptr<LineSearch> line_search_;
  std::string iteration_header_;
  //@}

}; // end Strategies

} // namespace FaRSA

#endif /* __FARSASTRATEGIES_HPP__ */
