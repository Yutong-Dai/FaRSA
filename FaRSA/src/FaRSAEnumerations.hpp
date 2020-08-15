// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __FARSAENUMERATIONS_HPP__
#define __FARSAENUMERATIONS_HPP__

namespace FaRSA
{

/** @name Enumerations */
//@{
/**
 * Report type enumerations
 */
enum ReportType
{
  R_SOLVER = 0,
  R_SUBSOLVER
};
/**
 * Report level enumerations
 */
enum ReportLevel
{
  R_BASIC = 0,
  R_PER_ITERATION,
  R_PER_INNER_ITERATION
};
//@}

} // namespace FaRSA

#endif /* __FARSAENUMERATIONS_HPP__ */
