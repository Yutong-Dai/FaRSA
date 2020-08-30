/*
 * File: testSparseMat.hpp
 * Author: Yutong Dai (rothdyt@gmail.com)
 * File Created: 2020-08-26 21:44
 * Last Modified: 2020-08-26 22:58
 * --------------------------------------------
 * Description:
 */

#ifndef __TESTSPARSEMAT_HPP__
#define __TESTSPARSEMAT_HPP__

#include <iostream>

#include "FaRSAReporter.hpp"
#include "FaRSASparseMat.hpp"

using namespace FaRSA;

int testSparseMatImplementation(int option) {
  int result = 0;
  Reporter reporter;
  if (option == 1) {
    std::shared_ptr<StreamReport> s(new StreamReport("s", R_SOLVER, R_BASIC));
    s->setStream(&std::cout);
    reporter.addReport(s);
  }
  std::valarray<double> val = {5.0, 8.0, 3.0, 6.0};
  std::vector<int> rowidx = {0, 0, 2, 3, 4};
  std::vector<int> colidx = {0, 1, 2, 1};
  csrMat m(4, 4, &val, &rowidx, &colidx);
  m.print(&reporter, "m");
  if (option == 1) {
    // Print final message
    if (result == 0) {
      reporter.printf(R_SOLVER, R_BASIC, "TEST WAS SUCCESSFUL.\n");
    } else {
      reporter.printf(R_SOLVER, R_BASIC, "TEST FAILED.\n");
    }
  }
  return result;
}

#endif /* __TESTSPARSEMAT_HPP__ */
