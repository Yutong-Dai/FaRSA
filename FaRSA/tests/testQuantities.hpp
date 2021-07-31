// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __TESTOPTIONS_HPP__
#define __TESTOPTIONS_HPP__

#include <iostream>

#include "FaRSAOptions.hpp"
#include "FaRSAReporter.hpp"
#include "FaRSAQuantities.hpp"

using namespace FaRSA;

// Implementation of test
int testQuantitiesImplementation(int option)
{

  // Initialize output
  int result = 0;

  // Declare reporter
  Reporter reporter;

  // Check option
  if (option == 1) {

    // Declare stream report
    std::shared_ptr<StreamReport> sr(new StreamReport("s", R_SOLVER, R_BASIC));

    // Set stream report to standard output
    sr->setStream(&std::cout);

    // Add stream report to reporter
    reporter.addReport(sr);

  } // end if

  // Declare option
  Options o;
  Quantities quantities;
  // Print empty option list
  reporter.printf(R_SOLVER, R_BASIC, "Printing options list... should be empty:\n");
  o.print(&reporter);
  // add default options
  quantities.addOptions(&o, &reporter);
  // Print option list
  reporter.printf(R_SOLVER, R_BASIC, "Printing options list... should not be empty:\n");
  o.print(&reporter);
  // load values from options
  quantities.getOptions(&o, &reporter);
  quantities.print(&reporter);


  // Check option
  if (option == 1) {
    // Print final message
    if (result == 0) {
      reporter.printf(R_SOLVER, R_BASIC, "TEST WAS SUCCESSFUL.\n");
    }
    else {
      reporter.printf(R_SOLVER, R_BASIC, "TEST FAILED.\n");
    }
  } // end if

  // Return
  return result;

} // end testOptionsImplementation

#endif /* __TESTOPTIONS_HPP__ */
