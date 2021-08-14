// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#include "FaRSASpacePartitionGroupL1PGBased.hpp"

#include <iostream>

#include "FaRSADeclarations.hpp"
#include "FaRSADefinitions.hpp"
namespace FaRSA
{
// Add options
void SpacePartitionGroupL1PGBased::addOptions(Options* options, const Reporter* reporter)
{
    // Add double options
    options->addDoubleOption(reporter, "SPGL1_complexity_p", 1.0, 0.0, FARSA_DOUBLE_INFINITY,
                             "A parameter controls the global iteation complexity result.\n"
                             "Default     : 1.0.");
    options->addDoubleOption(reporter, "SPGL1_kappa1", 1e-1, 0.0, FARSA_DOUBLE_INFINITY,
                             "A factor assocociated with the norm of gradient of the composite function restict to a "
                             "given group. Please refer to the FaRSA-Group paper for details.\n"
                             "Default     : 1e-1.");
    options->addDoubleOption(reporter, "SPGL1_kappa2", 1e-2, 0.0, FARSA_DOUBLE_INFINITY,
                             "A factor assocociated with the norm of gradient of the composite function restict to a "
                             "colleciton of groups. Please refer to the FaRSA-Group paper for details.\n"
                             "Default     : 1e-2.");
    options->addDoubleOption(
        reporter, "SPGL1_gamma", 1.0, 0.0, FARSA_DOUBLE_INFINITY,
        "A factor to scale the optimlaity measuer over the second order variables. If it is "
        "larger than 1.0, then algorithm tends to take more second order direction computation. Please refer "
        "to the FaRSA-Group paper for details.\n"
        "Default     : 1.0.");
}  // end addOptions

// Set options
void SpacePartitionGroupL1PGBased::getOptions(const Options* options, const Reporter* reporter)
{
    options->valueAsDouble(reporter, "SPGL1_complexity_p", p_);
    options->valueAsDouble(reporter, "SPGL1_kappa1", kappa1_);
    options->valueAsDouble(reporter, "SPGL1_kappa2", kappa2_);
    options->valueAsDouble(reporter, "SPGL1_gamma", gamma_);
}  // end getOptions

// Initialize
void SpacePartitionGroupL1PGBased::initialize(const Options* options, Quantities* quantities, const Reporter* reporter)
{
    std::shared_ptr<std::vector<int>> groups_full_indices(new std::vector<int>());
    for (int i = 0; i < quantities->numberOfGroups(); i++)
    {
        groups_full_indices->push_back(i);
    }
    groups_full_indices_ = groups_full_indices;

    std::shared_ptr<std::vector<int>> variables_full_indices(new std::vector<int>());
    for (int i = 0; i < quantities->numberOfVariables(); i++)
    {
        variables_full_indices->push_back(i);
    }
    variables_full_indices_ = variables_full_indices;
}

// Partition Space
void SpacePartitionGroupL1PGBased::partitionSpace(const Options* options, Quantities* quantities,
                                                  const Reporter* reporter, Strategies* strategies)
{
    setStatus(SP_UNSET);
    try
    {
        auto group_first_order = std::make_shared<std::vector<int>>();
        auto group_second_order = std::make_shared<std::vector<int>>();
        auto indicies_working = std::make_shared<std::vector<int>>();
        auto current_iterate = quantities->currentIterate();

        // get per group norm, gradient (over composite function) norm, proximal gradient update norm
        bool evaluation_success = current_iterate->evaluatePerGroupStatistics(*quantities);
        if (!evaluation_success)
        {
            THROW_EXCEPTION(SP_PGBASED_PARTITION_EXCEPTION,
                            "Space Partition failed when evaluatePerGroupStatistics was called.")
        }
        // compute working groups
        std::shared_ptr<std::vector<int>> groups_second_order(new std::vector<int>());

        // find all groups that are currently non-zero, predicted to be non-zero and sufficientlt far from zero
        // std::cout << "\nkappa1_: " << kappa1_ << "\n" << std::endl;
        for (int i = 0; i < quantities->numberOfGroups(); i++)
        {
            // std::cout << "per_group_2norm_: " << current_iterate->perGroup2Norm()->values()[i]
            //           << " per_group_gradient_: " << current_iterate->perGroupGradientAll2Norm()->values()[i]
            //           << " per_group_proximal_: " <<
            //           current_iterate->perGroupProximalGradientUpdate2Norm()->values()[i]
            //           << std::endl;
            if ((fabs(current_iterate->perGroup2Norm()->values()[i] - 0.0) > 1e-16) &&
                (fabs(current_iterate->perGroupProximalGradientUpdate2Norm()->values()[i] - 0.0) > 1e-16))
            {
                if ((current_iterate->perGroup2Norm()->values()[i] >=
                     kappa1_ * current_iterate->perGroupGradientAll2Norm()->values()[i]))
                {
                    groups_second_order->push_back(i);
                }
            }
        }

        // std::cout << "First screen:\n" << std::endl;
        // for (int i = 0; i < groups_second_order->size(); i++)
        // {
        //     std::cout << " groups_second_order: " << (*groups_second_order)[i] << std::endl;
        // }
        // std::cout << "=========" << std::endl;

        // secondary screeing remove groups that are potentially close
        std::shared_ptr<std::vector<int>> indicies_second_order(new std::vector<int>());
        if (groups_second_order->size() != 0)
        {
            double                     small_radius = 0.0;
            int                        groups_working_size = 0;
            std::vector<int>::iterator pos;
            for (auto i : *groups_second_order)
            {
                small_radius += pow(current_iterate->perGroupGradientAll2Norm()->values()[i], 2.0);
                groups_working_size += (*quantities->groups())[i][1] - (*quantities->groups())[i][0] + 1;
            }
            small_radius = pow(small_radius, 1.0 / p_);
            for (auto i : *groups_second_order)
            {
                int    group_size = (*quantities->groups())[i][1] - (*quantities->groups())[i][0] + 1;
                double scaling_factor = group_size / groups_working_size;
                if (current_iterate->perGroupGradientAll2Norm()->values()[i] <
                    fmin(kappa2_ * scaling_factor * small_radius, 1.0))
                {
                    pos = std::find(groups_second_order->begin(), groups_second_order->end(), i);
                    if (pos != groups_second_order->end())
                    {
                        groups_second_order->erase(pos);
                    }
                }
                else
                {
                    int start = (*quantities->groups())[i][0];
                    int end = (*quantities->groups())[i][1];
                    for (int j = start; j <= end; j++)
                    {
                        indicies_second_order->push_back(j);
                    }
                }
            }
        }
        // std::cout << "Secondary screen:\n" << std::endl;
        // for (int i = 0; i < groups_second_order->size(); i++)
        // {
        //     std::cout << " groups_second_order: " << (*groups_second_order)[i] << std::endl;
        // }
        // std::cout << "=========" << std::endl;

        // compute optimality measure
        auto   s = current_iterate->proximalGraidentStep();
        double s_norm2 = s->norm2();
        double chi_second_order = 0.0;
        double chi_first_order = 0.0;
        for (auto i : *indicies_second_order)
        {
            chi_second_order += pow(s->values()[i], 2.0);
        }
        chi_first_order = sqrt(pow(s_norm2, 2.0) - chi_second_order);
        chi_second_order = sqrt(chi_second_order);
        // std::cout << "\ns_norm2: " << s_norm2 << " chi_first_order: " << chi_first_order
        //           << " chi_second_order: " << chi_second_order << std::endl;
        // decide which direction should be computed
        if (chi_first_order <= gamma_ * chi_second_order)
        {
            strategies->directionComputationFirstOrder()->setPerformComputation(false);
            strategies->directionComputationSecondOrder()->setPerformComputation(true);
            // todo: allow to select a subset of groups
            auto groups_second_order_selected = groups_second_order;
            auto indicies_second_order_selected = indicies_second_order;
            quantities->setGroupsWorking(groups_second_order_selected);
            quantities->setIndiciesWorking(indicies_second_order_selected);
        }
        else
        {
            strategies->directionComputationFirstOrder()->setPerformComputation(true);
            strategies->directionComputationSecondOrder()->setPerformComputation(false);
            std::shared_ptr<std::vector<int>> groups_first_order(new std::vector<int>());
            std::shared_ptr<std::vector<int>> indicies_first_order(new std::vector<int>());

            // the following code will break if groups_second_order are not sorted; the same applies to
            // indicies_second_order compute groups_first_order
            std::set_difference(groups_full_indices_->begin(), groups_full_indices_->end(),
                                groups_second_order->begin(), groups_second_order->end(),
                                std::inserter(*groups_first_order, groups_first_order->begin()));
            // compute indicies_first_order
            std::set_difference(variables_full_indices_->begin(), variables_full_indices_->end(),
                                indicies_second_order->begin(), indicies_second_order->end(),
                                std::inserter(*indicies_first_order, indicies_first_order->begin()));
            quantities->setGroupsWorking(groups_first_order);
            quantities->setIndiciesWorking(indicies_first_order);
        }

        // Check for success
        THROW_EXCEPTION(SP_SUCCESS_EXCEPTION, "Space Partition successful.")
    }  // end try

    // catch exceptions
    catch (SP_SUCCESS_EXCEPTION& exec)
    {
        // Set status
        setStatus(SP_SUCCESS);
    }
    catch (SP_PGBASED_PARTITION_EXCEPTION& exec)
    {
        setStatus(SP_PG_BASED_PARTITION_FAILURE);
    }
}

}  // namespace FaRSA
