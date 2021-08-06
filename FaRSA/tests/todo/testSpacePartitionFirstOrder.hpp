// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __TESTOPTIONS_HPP__
#define __TESTOPTIONS_HPP__

#include <iostream>

#include "FaRSAOptions.hpp"
#include "FaRSAQuantities.hpp"
#include "FaRSAReporter.hpp"
#include "FaRSASpacePartitionFirstOrder.hpp"
using namespace FaRSA;

// Implementation of test
int testSpacePartitionFirstOrderImplementation(int option)
{
    // Initialize output
    int result = 0;

    // Declare reporter
    Reporter reporter;

    // Check option
    if (option == 1)
    {
        // Declare stream report
        std::shared_ptr<StreamReport> sr(
            new StreamReport("s", R_SOLVER, R_BASIC));

        // Set stream report to standard output
        sr->setStream(&std::cout);

        // Add stream report to reporter
        reporter.addReport(sr);

    }  // end if

    SpacePartitionFirstOrder partition;
    reporter.printf(R_SOLVER, R_BASIC, "Partition name: %s\n",
                    partition.name().c_str());

    Quantities quantities;
    Options    o;
    Strategies strategies;
    quantities.setNumberOfGroups(20);
    partition.partitionSpace(&o, &quantities, &reporter, &strategies);
    auto groups_first_order = quantities.groupsFirstOrder();
    reporter.printf(R_SOLVER, R_BASIC, "Groups in the groups_first_order:\n");
    for (int i = 0; i < groups_first_order->size(); i++)
    {
        reporter.printf(R_SOLVER, R_BASIC, "%2d ", (*groups_first_order)[i]);
        if ((i + 1) % 10 == 0)
        {
            reporter.printf(R_SOLVER, R_BASIC, "\n");
        }
    }
    for (int i = 0; i < groups_first_order->size(); i++)
    {
        if ((*groups_first_order)[i] != i)
        {
            result = 0;
            reporter.printf(R_SOLVER, R_BASIC,
                            "Expect to be %2d, but get %2d\n", i,
                            (*groups_first_order)[i]);
        }
    }

    // Check option
    if (option == 1)
    {
        // Print final message
        if (result == 0)
        {
            reporter.printf(R_SOLVER, R_BASIC, "TEST WAS SUCCESSFUL.\n");
        }
        else
        {
            reporter.printf(R_SOLVER, R_BASIC, "TEST FAILED.\n");
        }
    }  // end if

    // Return
    return result;

}  // end testOptionsImplementation

#endif /* __TESTOPTIONS_HPP__ */
