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
GroupL1::GroupL1(std::shared_ptr<std::vector<std::vector<int>>> group,
                 std::shared_ptr<std::vector<double>> weights, double penalty)
{
    name_ = "GroupL1";
    penalty_ = penalty;
    if ((*group).size() != (*weights).size())
    {
        THROW_EXCEPTION(FARSA_FUNCTION_INITIALIZATION_EXCEPTION,
                        "Group L1 initialization failed. Unmatched number of "
                        "groups inferred from the group and weights");
    }
    // declare a temp 2-d vector
    std::shared_ptr<std::vector<std::vector<int>>> g(
        new std::vector<std::vector<int>>);
    std::shared_ptr<std::vector<double>> w(new std::vector<double>);
    groups_ = g;
    weights_ = w;
    number_of_variables_ = 0;
    for (int i = 0; i < group->size(); i++)
    {
        std::vector<int> temp;
        for (int j = 0; j < ((*group)[i]).size(); j++)
        {
            temp.push_back((*group)[i][j]);
            idx_to_group_.push_back(i);
        }
        groups_->push_back(temp);
        // get the actual weights for each group
        weights_->push_back(((*weights)[i]) * penalty_);
        //
        number_of_variables_ += temp.size();
    }

    per_group_norm_ = std::vector<double>(groups_->size(), 0.0);
}

// evaluateObjective
bool GroupL1::evaluateObjective(const Vector& x, double& f)
{
    for (int i = 0; i < groups_->size(); i++)
    {
        double group_i_norm;
        group_i_norm = 0.0;
        for (int idx : (*groups_)[i])
        {
            group_i_norm += pow(x.values()[idx], 2);
        }
        group_i_norm = sqrt(group_i_norm);
        f = f + group_i_norm * (*weights_)[i];
        per_group_norm_[i] = group_i_norm;
    }

    return !isnan(f);
}  // end  evaluateObjective

// evaluateGradient
bool GroupL1::evaluateGradient(const Vector&           x,
                               const std::vector<int>& group_idx, Vector& g)
{
    int              count = 0;
    std::vector<int> cols;
    for (int idx : group_idx)
    {
        for (int element : (*groups_)[idx])
        {
            cols.push_back(element);
        }
    }
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