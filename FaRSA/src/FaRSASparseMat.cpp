/*
 * File: FaRSASparseMat.cpp
 * Author: Yutong Dai (rothdyt@gmail.com)
 * File Created: 2020-08-26 20:02
 * Last Modified: 2020-08-26 22:50
 * --------------------------------------------
 * Description:
 */

#include "FaRSASparseMat.hpp"

namespace FaRSA {
SparseMat::SparseMat(
    unsigned int nrows,
    unsigned int ncols,
    std::valarray<double>* val, /*how to make it accept int/double/bool*/
    std::vector<int>* rowidx,
    std::vector<int>* colidx)
    : nrows_(nrows),
      ncols_(ncols),
      val_(val),
      rowidx_(rowidx),
      colidx_(colidx) {
  int nnz = val->size();
  int totalElements = ncols_ * nrows_;
  sparsity_ = nnz / totalElements;
}
SparseMat::~SparseMat() {
  rowidx_->clear();
  colidx_->clear();
}

void csrMat::print(const Reporter* reporter,
                   std::string name) const {
  reporter->printf(R_SOLVER, R_BASIC, "The elements of %s:\n", name.c_str());
  for (int row = 0; row < nrows_; row++) {
    unsigned int row_start = (*rowidx_)[row];
    unsigned int row_end = (*rowidx_)[row + 1];
    if (row_end == row_start) {
      // the current row is all zero
      reporter->printf(R_SOLVER, R_BASIC, "  row [%d]:  all 0.\n", row);
    } else {
      reporter->printf(R_SOLVER, R_BASIC, "  row [%d]: ", row);
      reporter->printf(R_SOLVER, R_BASIC, "[ ");
      for (int col = 0; col < ncols_; col++) {
        if ((col >= row_start) && (col < row_end)) {
          if (col == ncols_ - 1) {
            reporter->printf(R_SOLVER, R_BASIC, "%3.2f ]", (*val_)[col]);
          } else {
            reporter->printf(R_SOLVER, R_BASIC, "%3.2f, ", (*val_)[col]);
          }
        } else {
          if (col == ncols_ - 1) {
            reporter->printf(R_SOLVER, R_BASIC, "0 ]");
          } else {
            reporter->printf(R_SOLVER, R_BASIC, "0, ");
          }
        }
      }
      reporter->printf(R_SOLVER, R_BASIC, "\n");
    }
    // reporter->printf(R_SOLVER, R_BASIC, "]\n");
  }
}
}  // namespace FaRSA
