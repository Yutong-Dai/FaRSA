// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#include "FaRSASpacePartitionFirstOrder.hpp"

#include "FaRSADeclarations.hpp"
#include "FaRSADefinitions.hpp"

namespace FaRSA
{
// Add options
void SpacePartitionFirstOrder::addOptions(Options*        options,
                                          const Reporter* reporter)
{
}  // end addOptions

// Set options
void SpacePartitionFirstOrder::getOptions(const Options*  options,
                                          const Reporter* reporter)
{
}  // end getOptions

// Initialize
void SpacePartitionFirstOrder::initialize(const Options*  options,
                                          Quantities*     quantities,
                                          const Reporter* reporter)
{
}

// Partition Space
void SpacePartitionFirstOrder::partitionSpace(const Options*  options,
                                              Quantities*     quantities,
                                              const Reporter* reporter,
                                              Strategies*     strategies)
{
    setStatus(SP_UNSET);
    try
    {
        // all variables should be seen as the first order variables
        auto group_first_order =
            std::make_shared<std::vector<int>>(quantities->numberOfGroups());
        auto group_second_order = std::make_shared<std::vector<int>>();

        for (int i = 0; i < group_first_order->size(); i++)
        {
            (*group_first_order)[i] = i;
        }
        quantities->setGroupsFirstOrder(group_first_order);
        quantities->setGroupsSecondOrder(group_second_order);
        if (quantities->groupsFirstOrder()->size() == 0 or
            quantities->groupsSecondOrder()->size() != 0)
        {
            THROW_EXCEPTION(
                SP_FIRST_ORDER_PARTITION_EXCEPTION,
                "Space Partition unsuccessful. The groups_first_order "
                "is empty or the groups_second_order is not empty.")
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
    catch (SP_FIRST_ORDER_PARTITION_EXCEPTION& exec)
    {
        setStatus(SP_FIRST_ORDER_PARTITION_FAILURE);
    }
}

}  // namespace FaRSA
