// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __FARSAPROBLEM_HPP__
#define __FARSAPROBLEM_HPP__

namespace FaRSA
{

/**
 * Problem class
 */
class Problem
{

public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   */
  Problem(){};
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  virtual ~Problem(){};

  /** @name Get methods */
  //@{
  /**
   * Returns number of variables
   * \param[out] n is the number of variables, an integer (return value)
   */
  virtual bool numberOfDataPoints(int& N) = 0;
  /**
   * Returns number of variables
   * \param[out] n is the number of variables, an integer (return value)
   */
  virtual bool numberOfVariables(int& n) = 0;
  /**
   * Returns initial point
   * \param[in] n is the number of variables, the size of "x", a constant integer
   * \param[out] x is the initial point/iterate, a double array (return value)
   */
  virtual bool initialPoint(int n,
                            double* x) = 0;
  //@}

  /** @name Evaluate methods */
  //@{
  /**
   * Evaluates objective
   * \param[in] n is the number of variables, the size of "x", a constant integer
   * \param[in] N is the number of data points, a constant integer
   * \param[in] x is a given point/iterate, a constant double array
   * \param[out] f is the objective value at "x", a double (return value)
   */
  virtual bool evaluateObjective(int n,
                                 int N,
                                 const double* x,
                                 double& f) = 0;
  /**
   * Evaluates gradient
   * \param[in] n is the number of variables, the size of "x", a constant integer
   * \param[in] N is the number of data points, a constant integer
   * \param[in] x is a given point/iterate, a constant double array
   * \param[out] g is the gradient value at "x", a double array (return value)
   */
  virtual bool evaluateGradient(int n,
                                int N,
                                const double* x,
                                double* g) = 0;
  /**
   * Evaluates Hessian-vector product
   * \param[in] n is the number of variables, the size of "x", a constant integer
   * \param[in] N is the number of data points, a constant integer
   * \param[in] x is a given point/iterate, a constant double array
   * \param[in] v is a given vector, a constant double array
   * \param[out] Hv is the Hessian value at "x" times "v", a double array (return value)
   */
  virtual bool evaluateHessianVectorProduct(int n,
                                            int N,
                                            const double* x,
                                            const double* v,
                                            double* Hv) = 0;
  //@}

  /** @name Finalize methods */
  //@{
  /**
   * Finalizes solution
   * \param[in] n is the number of variables, the size of "x", a constant integer
   * \param[in] N is the number of data points, a constant integer
   * \param[in] x is the final point/iterate, a constant double array
   * \param[in] f is the objective value at "x", a constant double
   * \param[in] g is the gradient value at "x", a constant double array
   */
  virtual bool finalizeSolution(int n,
                                int N,
                                const double* x,
                                double f,
                                const double* g) = 0;
  //@}

private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Copy constructor
   */
  Problem(const Problem&);
  /**
   * Overloaded equals operator
   */
  void operator=(const Problem&);
  //@}

}; // end Problem

} // namespace FaRSA

#endif /* __FARSAPROBLEM_HPP__ */
