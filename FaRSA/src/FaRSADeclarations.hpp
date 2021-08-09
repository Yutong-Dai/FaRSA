// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __FARSADECLARATIONS_HPP__
#define __FARSADECLARATIONS_HPP__

#include "FaRSAException.hpp"

namespace FaRSA
{
/** @name Exceptions */
//@{
/**
 * FaRSA exceptions
 */
DECLARE_EXCEPTION(FARSA_SUCCESS_EXCEPTION);

// For initialization
DECLARE_EXCEPTION(FARSA_INITIALIZATION_FAILURE_EXCEPTION);
/* this is used when instatiate a function object */
/* it will be caught outside of the solver class */
DECLARE_EXCEPTION(FARSA_FUNCTION_INITIALIZATION_EXCEPTION);

// for termination - based on limits
DECLARE_EXCEPTION(FARSA_CPU_TIME_LIMIT_EXCEPTION);
DECLARE_EXCEPTION(FARSA_ITERATE_NORM_LIMIT_EXCEPTION);
DECLARE_EXCEPTION(FARSA_ITERATION_LIMIT_EXCEPTION);
DECLARE_EXCEPTION(FARSA_FUNCTION_EVALUATION_LIMIT_EXCEPTION);
DECLARE_EXCEPTION(FARSA_GRADIENT_EVALUATION_LIMIT_EXCEPTION);
DECLARE_EXCEPTION(FARSA_HESSIAN_VECTOR_PRODUCT_EVALUATION_LIMIT_EXCEPTION);

// for evaluation/computation failure
DECLARE_EXCEPTION(FARSA_FUNCTION_EVALUATION_FAILURE_EXCEPTION);
DECLARE_EXCEPTION(FARSA_GRADIENT_EVALUATION_FAILURE_EXCEPTION);
DECLARE_EXCEPTION(FARSA_PROXIMAL_GRADIENT_COMPUTATION_FAILURE_EXCEPTION);
DECLARE_EXCEPTION(FARSA_HESSIAN_VECTOR_PRODUCT_EVALUATION_FAILURE_EXCEPTION);

// for assert
DECLARE_EXCEPTION(FARSA_FUNCTION_EVALUATION_ASSERT_EXCEPTION);
DECLARE_EXCEPTION(FARSA_GRADIENT_EVALUATION_ASSERT_EXCEPTION);
DECLARE_EXCEPTION(FARSA_PROXIMAL_GRADIENT_COMPUTATION_ASSERT_EXCEPTION);
DECLARE_EXCEPTION(FARSA_HESSIAN_VECTOR_PRODUCT_EVALUATION_ASSERT_EXCEPTION);

DECLARE_EXCEPTION(FARSA_MATRIX_EXCEPTION);
DECLARE_EXCEPTION(FARSA_MATRIX_ASSERT_EXCEPTION);
DECLARE_EXCEPTION(FARSA_VECTOR_EXCEPTION);
DECLARE_EXCEPTION(FARSA_VECTOR_ASSERT_EXCEPTION);

// For strategies
DECLARE_EXCEPTION(FARSA_SPACE_PARTITION_FAILURE_EXCEPTION);
DECLARE_EXCEPTION(FARSA_DIRECTION_COMPUTATION_FAILURE_EXCEPTION);
DECLARE_EXCEPTION(FARSA_LINE_SEARCH_FAILURE_EXCEPTION);
DECLARE_EXCEPTION(FARSA_PARAMETER_UPDATE_FAILURE_EXCEPTION);
/**
 * Direction computation exceptions
 */
DECLARE_EXCEPTION(SP_SUCCESS_EXCEPTION);
DECLARE_EXCEPTION(SP_FIRST_ORDER_PARTITION_EXCEPTION);
DECLARE_EXCEPTION(SP_PGBASED_PARTITION_EXCEPTION);
/**
 * Direction computation exceptions
 */
DECLARE_EXCEPTION(DC_SUCCESS_EXCEPTION);
DECLARE_EXCEPTION(DC_EVALUATION_FAILURE_EXCEPTION);
DECLARE_EXCEPTION(DC_ITERATION_LIMIT_EXCEPTION);
/**
 * Line search exceptions
 */
DECLARE_EXCEPTION(LS_SUCCESS_EXCEPTION);
DECLARE_EXCEPTION(LS_EVALUATION_FAILURE_EXCEPTION);
DECLARE_EXCEPTION(LS_STEPSIZE_TOO_SMALL_EXCEPTION);
DECLARE_EXCEPTION(LS_ITERATION_LIMIT_EXCEPTION);
DECLARE_EXCEPTION(LS_MAXBACKTRACK_LIMIT_EXCEPTION);

/**
 * Parameter Update exceptions
 */
DECLARE_EXCEPTION(PU_SUCCESS_EXCEPTION);
DECLARE_EXCEPTION(PU_EVALUATION_FAILURE_EXCEPTION);
//@}

}  // namespace FaRSA

#endif /* __FARSADECLARATIONS_HPP__ */
