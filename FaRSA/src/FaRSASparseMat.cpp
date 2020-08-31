/*
 * File: FaRSASparseMat.cpp
 * Author: Yutong Dai (rothdyt@gmail.com)
 * File Created: 2020-08-26 20:02
 * Last Modified: 2020-08-31 11:56
 * --------------------------------------------
 * Description:
 */

#include "FaRSASparseMat.hpp"

#include <cassert>
#include <cstdio>
namespace FaRSA {
SparseMat::SparseMat(
    unsigned int nrows,
    unsigned int ncols,
    std::valarray<double> val, /*how to make it accept int/double/bool*/
    std::vector<int> rowidx,
    std::vector<int> colidx)
    : nrows_(nrows),
      ncols_(ncols),
      val_(val),
      rowidx_(rowidx),
      colidx_(colidx) {
  int nnz = val_.size();
  int totalElements = ncols_ * nrows_;
  sparsity_ = (float)nnz / (float)totalElements;  // be aware of integer divison
}
SparseMat::~SparseMat() {
  rowidx_.clear();
  colidx_.clear();
}

SparseMat SparseMat::operator*(double scalar) {
  std::valarray<double> newval_ = scalar * val_;
  return SparseMat(nrows_, ncols_, newval_, rowidx_, colidx_);
}

void csrMat::print(const Reporter* reporter,
                   std::string name) const {
  reporter->printf(R_SOLVER, R_BASIC, "The elements of %s:\n", name.c_str());
  for (int row = 0; row < nrows_; row++) {
    // printf("Working at row:[%d]\n", row);
    unsigned int row_start = rowidx_[row];
    unsigned int row_end = rowidx_[row + 1];
    // printf(" row_start:[%d]| row_end:[%d] \n", row_start, row_end);
    if (row_end == row_start) {
      // the current row is all zero
      reporter->printf(R_SOLVER, R_BASIC, "  row [%d]:  all 0.\n", row);
    } else {
      std::vector<float> valSub;
      std::vector<float> colSub;
      for (int i = row_start; i < row_end; i++) {
        valSub.push_back(val_[i]);
        colSub.push_back(colidx_[i]);
      }
      reporter->printf(R_SOLVER, R_BASIC, "  row [%d]: ", row);
      reporter->printf(R_SOLVER, R_BASIC, "[ ");
      int matched = 0;
      for (int col = 0; col < ncols_; col++) {
        bool condition = std::find(colSub.begin(), colSub.end(), col) != colSub.end();
        if (condition) {
          if (col == ncols_ - 1) {
            reporter->printf(R_SOLVER, R_BASIC, "%3.2f ]", valSub[matched]);
          } else {
            reporter->printf(R_SOLVER, R_BASIC, "%3.2f, ", valSub[matched]);
          }
          matched += 1;
        } else {
          if (col == ncols_ - 1) {
            reporter->printf(R_SOLVER, R_BASIC, "0.00 ]");
          } else {
            reporter->printf(R_SOLVER, R_BASIC, "0.00, ");
          }
        }
      }
      reporter->printf(R_SOLVER, R_BASIC, "\n");
    }
  }
}

std::shared_ptr<Vector> csrMat::dot(const Vector& v) {
  assert(v.length() == ncols_);
  std::shared_ptr<Vector>
      ans(new Vector(v.length()));
  for (int row = 0; row < nrows_; row++) {
    // printf("Working at row:[%d]\n", row);
    unsigned int row_start = rowidx_[row];
    unsigned int row_end = rowidx_[row + 1];
    // printf(" row_start:[%d]| row_end:[%d] \n", row_start, row_end);
    if (row_end == row_start) {
      ans->set(row, 0.0);
    } else {
      double temp = 0;
      for (int i = row_start; i < row_end; i++) {
        int colId = colidx_[i];
        temp += val_[i] * v.values()[colId];
        // printf("i:[%d] | colId:[%d] | matrix value:[%f] | vector value:[%f]\n", i, colId, val_[i], v.values()[colId]);
      }
      ans->set(row, temp);
    }
  }
  return ans;
}

// Vector csrMat::dot(const Vector& v) {
//   Vector ans(v.length());
//   for (int row = 0; row < nrows_; row++) {
//     unsigned int row_start = rowidx_[row];
//     unsigned int row_end = rowidx_[row + 1];
//     if (row_end == row_start) {
//       ans.set(row, 0.0);
//     } else {
//       double temp = 0;
//       for (int i = row_start; i < row_end; i++) {
//         int colId = colidx_[i];
//         temp += val_[i] * v.values()[colId];
//       }
//       ans.set(row, temp);
//     }
//   }
//   return ans;
// }
}  // namespace FaRSA
