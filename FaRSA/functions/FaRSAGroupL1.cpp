// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#include "FaRSAGroupL1.hpp"

// #include <iostream>

#include "FaRSADeclarations.hpp"
#include "FaRSADefinitions.hpp"
// constructor
GroupL1::GroupL1(std::vector<std::vector<int>>& groups,
                 std::vector<double>& weights, double penalty)
{
    name_ = "GroupL1";
    penalty_ = penalty;
    groups_ = groups;
    weights_ = weights;
    // get the actual weights for each group
    number_of_variables_ = 0;
    for (int i = 0; i < weights.size(); i++)
    {
        weights_[i] = weights_[i] * penalty_;
        number_of_variables_ =
            number_of_variables_ + (groups_[i][1] - groups_[i][0] + 1);
    }
    idx_to_group_.resize(number_of_variables_);
    for (int i = 0; i < groups_.size(); i++)
    {
        int start = groups_[i][0];
        int end = groups_[i][1];
        for (int j = start; j <= end; j++)
        {
            idx_to_group_[j] = i;
        }
    }
    per_group_norm_ = std::vector<double>(groups_.size(), 0.0);
}

// evaluateObjective
bool GroupL1::evaluateObjective(const Vector& x, double& f)
{
    for (int i = 0; i < groups_.size(); i++)
    {
        double group_i_norm;
        group_i_norm = 0.0;
        int start = groups_[i][0];
        int end = groups_[i][1];
        for (int j = start; j <= end; j++)
        {
            group_i_norm = group_i_norm + pow(x.values()[j], 2);
        }
        group_i_norm = sqrt(group_i_norm);
        f = f + group_i_norm * weights_[i];
        per_group_norm_[i] = group_i_norm;
    }

    return !isnan(f);
}  // end  evaluateObjective

// evaluateGradient
bool GroupL1::evaluateGradient(const Vector& x, const std::vector<int>& cols,
                               Vector& g)
{
    int group_idx;
    int count = 0;
    for (int i : cols)
    {
        group_idx = idx_to_group_[i];
        g.valuesModifiable()[count] =
            (weights_[group_idx] / per_group_norm_[group_idx]) * x.values()[i];
        if (isnan(g.values()[i]))
        {
            return false;
        }
        count = count + 1;
    }
    return true;
}  // end  evaluateGradient

// evaluateHessianVectorProduct
/*
    weights_i( V_{Gi}/||X_{Gi}|| - (X_{Gi}.dot(V_{Gi})/||X_{Gi}||^3 * X_{Gi}))
*/
bool GroupL1::evaluateHessianVectorProduct(const Vector&           x,
                                           const std::vector<int>& cols,
                                           const Vector& v, Vector& Hv)
{
    int                   group_idx;
    std::map<int, double> xgi_vgi_innerproduct;
    int                   count = 0;
    for (int i = 0; i < cols.size(); i++)
    {
        group_idx = idx_to_group_[cols[i]];
        std::map<int, double>::iterator iter =
            xgi_vgi_innerproduct.find(group_idx);
        if (iter == xgi_vgi_innerproduct.end())
        {
            // <XGi,VGi> is not computed
            double inner_prod_Gi = 0.0;
            int    start = groups_[group_idx][0];
            int    end = groups_[group_idx][1];
            for (int j = start; j <= end; j++)
            {
                inner_prod_Gi =
                    inner_prod_Gi + x.values()[j] * v.values()[count];
                count += 1;
            }
            xgi_vgi_innerproduct[group_idx] = inner_prod_Gi;
        }

        Hv.valuesModifiable()[i] =
            weights_[group_idx] *
            (v.values()[i] / per_group_norm_[group_idx] -
             xgi_vgi_innerproduct[group_idx] * x.values()[cols[i]] /
                 pow(per_group_norm_[group_idx], 3));

        if (isnan(Hv.values()[i]))
        {
            return false;
        }
    }
    return true;
}
// end  evaluateHessianVectorProduct

// evaluateProximalGradient
bool GroupL1::evaluateProximalGradient(const Vector& x, const Vector& gradfx,
                                       double stepsize, Vector& proxgrad)
{
    assert(proxgrad.length() == x.length());
    assert(proxgrad.length() == gradfx.length());
    Vector gradient_step(x.length());
    gradient_step.linearCombination(1.0, x, -stepsize, gradfx);
    int group_idx;
    per_group_gradient_step_norm_.clear();
    for (int i = 0; i < gradient_step.length(); i++)
    {
        group_idx = idx_to_group_[i];
        std::map<int, double>::iterator iter =
            per_group_gradient_step_norm_.find(group_idx);
        if (iter == per_group_gradient_step_norm_.end())
        {
            //  ||X_{Gi}-stepsize * gradfx_{Gi}|| is not computed
            double gradient_step_Gi_norm = 0.0;
            int    start = groups_[group_idx][0];
            int    end = groups_[group_idx][1];
            for (int j = start; j <= end; j++)
            {
                gradient_step_Gi_norm =
                    gradient_step_Gi_norm + pow(gradient_step.values()[j], 2);
            }

            per_group_gradient_step_norm_[group_idx] =
                sqrt(gradient_step_Gi_norm);
        }
        proxgrad.valuesModifiable()[i] =
            fmax(0, 1 - weights_[group_idx] * stepsize /
                            per_group_gradient_step_norm_[group_idx]) *
            gradient_step.values()[i];
        if (isnan(proxgrad.values()[i]))
        {
            return false;
        }
    }
    return true;
}