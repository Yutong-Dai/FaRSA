/*
 * File: testSparseMat.hpp
 * Author: Yutong Dai (rothdyt@gmail.com)
 * File Created: 2020-08-26 21:44
 * Last Modified: 2020-08-30 19:46
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
  std::vector<int> colidx = {1, 3, 0, 2};
  csrMat m(4, 4, val, rowidx, colidx);
  reporter.printf(R_SOLVER, R_BASIC, "The sparsity of m is %f\n", m.sparsity());
  m.print(&reporter, "m");

  csrMat mnew = m * 3.0;  // be careful of type conversion
  reporter.printf(R_SOLVER, R_BASIC, "The sparsity of mnew is %f\n", mnew.sparsity());
  mnew.print(&reporter, "m new");

  Vector v(4);
  v.set(0, 1);
  v.set(1, 2);
  v.set(2, 3);
  v.set(3, 4);
  v.print(&reporter, "Vector v:");
  std::shared_ptr<Vector> ans = m.dot(v);

  for (int i = 0; i < 4; i++) {
    reporter.printf(R_SOLVER, R_BASIC, "The [%d]-th element is %f\n", i, ans->values()[i]);
  }

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
