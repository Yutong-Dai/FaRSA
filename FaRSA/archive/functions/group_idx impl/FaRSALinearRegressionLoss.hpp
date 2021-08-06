// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

// Description : Implementation for FaRSA::LinearRegressionLoss
//                 f(x;A,b) = 1/(2m) * ||Ax-b||^2, A\in R^{m,n}, x\in R^{n}.

#ifndef __LINEARREGRESSIONLOSS_HPP__
#define __LINEARREGRESSIONLOSS_HPP__

#include <vector>

#include "FaRSADeclarations.hpp"
#include "FaRSADefinitions.hpp"
#include "FaRSAFunctionSmooth.hpp"
using namespace FaRSA;

class LinearRegressionLoss : public FunctionSmooth
{
   public:
    /** @name Constructors */
    //@{
    /**
     * Constructor
     */
    LinearRegressionLoss(char*            dataset_file_path,
                         SparseFormatType sparse_format, char* label_file_path,
                         std::string dataset_name,
                         std::shared_ptr<std::vector<std::vector<int>>> group)
    {
        data_matrix_.setFromFile(dataset_file_path, sparse_format);
        data_label_.setFromFile(label_file_path);
        number_of_data_points_ = data_matrix_.numberOfRows();
        number_of_variables_ = data_matrix_.numberOfColumns();
        residual_.setLength(data_matrix_.numberOfRows());
        dataset_name_ = dataset_name;
        name_ = "LinearRegressionLoss";
        groups_ = group;
    };
    //@}

    /** @name Destructor */
    //@{
    /**
     * Destructor
     */
    ~LinearRegressionLoss(){};

    /** @name Get methods */
    //@{
    //@}

    /** @name Evaluate methods */
    //@{
    /**
     * Evaluates objective
     * \param[in] x is a given point/iterate, a constant double array
     * \param[out] f is the objective value at "x", a double (return value)
     * \return indicator of success (true) or failure (false)
     */
    bool evaluateObjective(const Vector& x, double& f);
    /**
     * Evaluates gradient in subspace defined by variables
     * in subgroups.
     * \param[in] x is a given point/iterate, a constant double array
     * \param[in] group_idx is a vector of column indices; if empty, then
     * calculate full gradient.
     * \param[in] v is a given vector, a constant double
     * array \param[out] Hv is the product of the Hessian and "v", a double
     * array (return value) \return indicator of success (true) or failure
     * (false)
     */
    bool evaluateGradient(const Vector& x, const std::vector<int>& group_idx,
                          Vector& g);
    /**
     * Evaluates Hessian-vector product in subspace defined by variables
     * in cols.
     * \param[in] x is a given point/iterate, a constant double array
     * \param[in] group_idx is a vector of selected group indices
     * \param[in] v is a given vector, a constant double array
     * \param[out] Hv is the Hessian value at "x" times "v", a double array
     * (return value)
     */
    bool evaluateHessianVectorProduct(const Vector&           x,
                                      const std::vector<int>& group_idx,
                                      const Vector& v, Vector& Hv);
    //@}

   private:
    /** @name Default compiler generated methods
     * (Hidden to avoid implicit creation/calling.)
     */
    //@{
    /**
     * Copy constructor
     */
    LinearRegressionLoss(const LinearRegressionLoss&);
    /**
     * Overloaded equals operator
     */
    void operator=(const LinearRegressionLoss&);
    //@}
    /**
     * Private members
     */
    //@{
    Vector residual_; /* Ax-b */
    //@}

};  // end LinearRegressionLoss

#endif /* __LINEARREGRESSIONLOSS_HPP__ */
