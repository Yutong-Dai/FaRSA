// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#include <cstdio>

#include "testGroupL1.hpp"
#include "testLinearRegressionLoss.hpp"
#include "testPoint.hpp"
// Main function
int main()
{
    // Initialize result
    int result = 0;

    // Run tests
    printf("testing GroupL1...................... ");
    if (!testGroupL1Implementation(0))
    {
        printf("success.\n");
    }
    else
    {
        result = 1;
        printf("failure! (run testGroupL1 for details)\n");
    }
    printf("testing LinearRegressionLoss......... ");
    if (!testLinearRegressionLossImplementation(0))
    {
        printf("success.\n");
    }
    else
    {
        result = 1;
        printf("failure! (run testLinearRegressionLoss for details)\n");
    }
    printf("testing Point........................ ");
    if (!testPointImplementation(0))
    {
        printf("success.\n");
    }
    else
    {
        result = 1;
        printf("failure! (run testPoint for details)\n");
    }

    // Return
    return 0;

}  // end main
