// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __TESTMATRIX_HPP__
#define __TESTMATRIX_HPP__

#include <iostream>
#include <vector>

#include "FaRSAEnumerations.hpp"
#include "FaRSAReporter.hpp"
#include "FaRSAMatrix.hpp"
#include "FaRSAVector.hpp"

using namespace FaRSA;

// Implementation of test
int testMatrixImplementation(int option)
{

  // Initialize output
  int result = 0;

  // Declare reporter
  Reporter reporter;

  // Check option
  if (option == 1) {

    // Declare stream report
    std::shared_ptr<StreamReport> s(new StreamReport("s", R_SOLVER, R_BASIC));

    // Set stream report to standard output
    s->setStream(&std::cout);

    // Add stream report to reporter
    reporter.addReport(s);

  } // end if

  // Declare matrix
  Matrix A;

  // Set file name
  char* file_name = (char*)"matrix.txt";

  // Read from file
  A.setFromFile(file_name, M_COORDINATE_LIST);

  // Print matrix
  A.print(&reporter,"Testing read from file:");

  // Create vectors for matrix-vector product
  Vector x(6);
  Vector b(3);

  // Set elements of vector
  x.valuesModifiable()[0] =  66.6;
  x.valuesModifiable()[1] = -55.5;
  x.valuesModifiable()[2] =  44.4;
  x.valuesModifiable()[3] = -33.3;
  x.valuesModifiable()[4] =  22.2;
  x.valuesModifiable()[5] = -11.1;

  // Compute matrix-vector product
  A.matrixVectorProduct(x,b);

  // Check values
  for (int i = 0; i < 3; i++) {
    if (b.values()[i] < -1e-12 || b.values()[i] > 1e-12) {
      result = 1;
    }
  } // end for

  // Print product
  b.print(&reporter,"Testing matrix-vector product:");

  // Create vectors for matrix-transpose-vector product
  Vector y(3);
  Vector c(6);

  // Set elements of vector
  y.valuesModifiable()[0] =  123.4;
  y.valuesModifiable()[1] = -432.1;
  y.valuesModifiable()[2] =  121.2;

  // Compute matrix-transpose-vector product
  A.matrixTransposeVectorProduct(y,c);

  // Check values
  if (c.values()[0] < 1.357400000000000e+02 - 1e-12 || c.values()[0] > 1.357400000000000e+02 + 1e-12) {
    result = 1;
  }
  if (c.values()[1] < -9.506200000000001e+02 - 1e-12 || c.values()[1] > -9.506200000000001e+02 + 1e-12) {
    result = 1;
  }
  if (c.values()[2] < 3.999600000000000e+02 - 1e-12 || c.values()[2] > 3.999600000000000e+02 + 1e-12) {
    result = 1;
  }
  if (c.values()[3] < 5.332800000000001e+02 - 1e-12 || c.values()[3] > 5.332800000000001e+02 + 1e-12) {
    result = 1;
  }
  if (c.values()[4] < -2.376550000000000e+03 - 1e-12 || c.values()[4] > -2.376550000000000e+03 + 1e-12) {
    result = 1;
  }
  if (c.values()[5] < 8.144399999999999e+02 - 1e-12 || c.values()[5] > 8.144399999999999e+02 + 1e-12) {
    result = 1;
  }

  // Print product
  c.print(&reporter,"Testing matrix-transpose-vector product:");

  // check col method
  std::vector<int>  col_idx = {0,2};
  Matrix Asub(A.numberOfRows(), col_idx.size(), col_idx.size() * A.numberOfRows());
  A.col(col_idx, Asub);
  Asub.print(&reporter,"Submatrix of A:");
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
  reporter.printf(R_SOLVER, R_BASIC, "End of the file about to call destructors.");
  // Return
  return result;

} // end testMatrixImplementation

#endif /* __TESTMATRIX_HPP__ */
