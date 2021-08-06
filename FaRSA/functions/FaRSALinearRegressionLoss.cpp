// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

// Description : Implementation for FaRSA::LinearRegressionLoss
//                 f(x;A,b) = 1/(2m) * ||Ax-b||^2, A\in R^{m,n}, x\in R^{n}.

#include "FaRSALinearRegressionLoss.hpp"

#include <iostream>

bool LinearRegressionLoss::evaluateObjective(const Vector& x, double& f)
{
    data_matrix_.matrixVectorProduct(x, residual_);
    residual_.addScaledVector(-1, data_label_);
    f = residual_.innerProduct(residual_) / (2 * number_of_data_points_);
    return !isnan(f);
}  // end  evaluateObjective

bool LinearRegressionLoss::evaluateGradient(const Vector&           x,
                                            const std::vector<int>& indicies,
                                            Vector&                 g)
{
    data_matrix_.matrixTransposeVectorProduct(residual_, g);
    g.scale(1.0 / number_of_data_points_);
    for (int i = 0; i < g.length(); i++)
    {
        if (isnan(g.values()[i]))
        {
            return false;
        }
    }
    return true;
}  // end  evaluateGradient

bool LinearRegressionLoss::evaluateHessianVectorProduct(
    const Vector& x, const std::vector<int>& indicies, const Vector& v,
    Vector& Hv)
{
    Matrix submatrix(data_matrix_.numberOfRows(), indicies.size(),
                     indicies.size() * data_matrix_.numberOfRows());
    data_matrix_.col(indicies, submatrix);
    Vector temp(submatrix.numberOfRows());
    submatrix.matrixVectorProduct(v, temp);
    submatrix.matrixTransposeVectorProduct(temp, Hv);
    Hv.scale(1.0 / number_of_data_points_);
    for (int i = 0; i < Hv.length(); i++)
    {
        if (isnan(Hv.values()[i]))
        {
            return false;
        }
    }
    return true;
}
// end  evaluateHessianVectorProduct