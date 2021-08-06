// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#include "FaRSAGroupL1.hpp"

#include "FaRSADeclarations.hpp"
#include "FaRSADefinitions.hpp"

// constructor
GroupL1::GroupL1(std::vector<std::vector<int>>& groups,
                 std::vector<float>& weights, float penalty)
{
    penalty_ = penalty;
    groups_ = groups;
    weights_ = weights;
    // get the actual weights for each group
    for (int i = 0; i < weights.size(); i++)
    {
        weights_[i] = weights_[i] * penalty_;
    }
    per_group_norm_ = std::vector<float>(groups_.size(), 0.0);
    for (int i = 0; i < groups_.size(); i++)
    {
        for (int j = 0; j < groups_[i].size(); j++)
        {
            idx_to_group_[groups_[i][j]] = i;
        }
    }
}

// evaluateObjective
bool GroupL1::evaluateObjective(const Vector& x, double& f)
{
    for (int i = 0; i < groups_.size(); i++)
    {
        double group_i_norm;
        group_i_norm = 0.0;
        for (int j = 0; j < groups_[i].size(); j++)
        {
            group_i_norm = group_i_norm + pow(x.values()[groups_[i][j]], 2);
        }
        f = f + sqrt(group_i_norm) * weights_[i];
        per_group_norm_[i] = group_i_norm;
    }

    return !isnan(f);
}  // end  evaluateObjective

// evaluateGradient
bool GroupL1::evaluateGradient(const Vector& x, const std::vector<int>& cols,
                               Vector& g)
{
    int group_idx;
    for (int i = 0; i < cols.size(); i++)
    {
        group_idx = idx_to_group_[cols[i]];
        g.valuesModifiable()[cols[i]] =
            (weights_[group_idx] / per_group_norm_[group_idx]) *
            x.values()[cols[i]];
        if (isnan(g.values()[cols[i]]))
        {
            return false;
        }
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
    int                  group_idx;
    std::map<int, float> xgi_vgi_innerproduct;
    for (int i = 0; i < cols.size(); i++)
    {
        group_idx = idx_to_group_[cols[i]];
        std::map<int, float>::iterator iter =
            xgi_vgi_innerproduct.find(group_idx);
        if (iter == xgi_vgi_innerproduct.end())
        {
            // <XGi,VGi> is not computed
            float inner_prod_Gi = 0.0;
            for (int idx : groups_[group_idx])
            {
                inner_prod_Gi =
                    inner_prod_Gi + x.values()[idx] * v.values()[idx];
            }
            xgi_vgi_innerproduct[group_idx] = inner_prod_Gi;
            Hv.valuesModifiable()[cols[i]] =
                weights_[group_idx] *
                (v.values()[cols[i]] / per_group_norm_[group_idx] -
                 inner_prod_Gi * x.values()[cols[i]] /
                     pow(per_group_norm_[group_idx], 3));
        }
        else
        {
            Hv.valuesModifiable()[cols[i]] =
                weights_[group_idx] *
                (v.values()[cols[i]] / per_group_norm_[group_idx] -
                 xgi_vgi_innerproduct[group_idx] * x.values()[cols[i]] /
                     pow(per_group_norm_[group_idx], 3));
        }
        if (isnan(Hv.values()[cols[i]]))
        {
            return false;
        }
    }
    return true;
}
// end  evaluateHessianVectorProduct

// evaluateProximalGradient
bool GroupL1::evaluateProximalGradient(const Vector& x, const Vector& gradfx,
                                       float stepsize, Vector& proxgrad)
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
        std::map<int, float>::iterator iter =
            per_group_gradient_step_norm_.find(group_idx);
        if (iter == per_group_gradient_step_norm_.end())
        {
            //  ||X_{Gi}-stepsize * gradfx_{Gi}|| is not computed
            float gradient_step_Gi_norm = 0.0;
            for (int idx : groups_[group_idx])
            {
                gradient_step_Gi_norm =
                    gradient_step_Gi_norm + pow(gradient_step.values()[idx], 2);
            }
            per_group_gradient_step_norm_[group_idx] =
                sqrt(gradient_step_Gi_norm);
            proxgrad.valuesModifiable()[i] =
                fmax(0, 1 - weights_[group_idx] * stepsize /
                                per_group_gradient_step_norm_[group_idx]) *
                gradient_step.values()[i];
        }
        else
        {
            proxgrad.valuesModifiable()[i] =
                fmax(0, 1 - weights_[group_idx] * stepsize /
                                per_group_gradient_step_norm_[group_idx]) *
                gradient_step.values()[i];
        }
        if (isnan(proxgrad.values()[i]))
        {
            return false;
        }
    }
    return true;
}