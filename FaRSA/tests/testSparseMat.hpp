/*
 * File: testSparseMat.hpp
 * Author: Yutong Dai (rothdyt@gmail.com)
 * File Created: 2020-08-26 21:44
 * Last Modified: 2020-09-28 01:31
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
  double val[] = {5.0, 8.0, 3.0, 6.0, 7.0};
  int rowPtr[] = {0, 0, 2, 3, 5};
  int colidx[] = {1, 3, 0, 2, 4};
  csrMat m(4, 5, val, rowPtr, colidx);
  m.print(&reporter, "m");

  Vector v(5);
  v.set(0, 1);
  v.set(1, 2);
  v.set(2, 3);
  v.set(3, 4);
  v.set(4, 5);
  v.print(&reporter, "Vector v:");
  std::shared_ptr<Vector> ans = m.dot(v);

  for (int i = 0; i < 4; i++) {
    reporter.printf(R_SOLVER, R_BASIC, "The [%d]-th element is %f\n", i, ans->values()[i]);
  }

  std::shared_ptr<SparseMat> mT = m.transpose();
  mT->print(&reporter, "m.T");
  Vector u(4);
  u.set(0, 1);
  u.set(1, 2);
  u.set(2, 3);
  u.set(3, 4);
  u.print(&reporter, "Vector u:");
  std::shared_ptr<Vector> ansT = mT->dot(u);
  for (int i = 0; i < 5; i++) {
    reporter.printf(R_SOLVER, R_BASIC, "The [%d]-th element is %f\n", i, ansT->values()[i]);
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
