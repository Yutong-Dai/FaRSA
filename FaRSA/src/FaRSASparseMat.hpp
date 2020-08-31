/*
 * File: FaRSASparseMat.h
 * Author: Yutong Dai (rothdyt@gmail.com)
 * File Created: 2020-08-26 15:38
 * Last Modified: 2020-08-30 18:57
 * --------------------------------------------
 * Description:
 */

#ifndef __FARSASPARSEMAT_HPP__
#define __FARSASPARSEMAT_HPP__

#include <memory>
#include <string>
#include <valarray>

#include "FaRSAReporter.hpp"
#include "FaRSAVector.hpp"
namespace FaRSA {

/**
 * Forward declarations
 */
class Reporter;

/**
 * Abstract Sparse Matrix class
 */
class SparseMat {
 public:
  /**@name Constructor*/
  //@{
  /**
   * Constructor; no argument
   */
  SparseMat()
      : nrows_(0),
        ncols_(0),
        val_({}),
        rowidx_(0),
        colidx_(0),
        sparsity_(-1.0){};
  /**
   * Constructor with given val, rowidx, colidx
   * \param[in] val is a pointer to an array with all nonzero elements in the SparseMat
   * \param[in] rowidx is a pointer to an array with row indices of nonzero elements
   * \param[in] colidx is a pointer to an array with column indices of nonzero elements
   */
  SparseMat(
      unsigned int nrows,
      unsigned int ncols,
      std::valarray<double> val,
      std::vector<int> rowidx,
      std::vector<int> colidx);
  SparseMat(const SparseMat& other)
      : nrows_(other.nrows_),
        ncols_(other.ncols_),
        val_(other.val_),
        rowidx_(other.rowidx_),
        colidx_(other.colidx_),
        sparsity_(other.sparsity_){};
  //@}
  /**@name Destructor*/
  //@{
  /** Destructor; val_, rowidx_, colidx_ delete*/
  ~SparseMat();
  //@}
  //*@name Print  method*/
  //@{
  /**
   * virtual method for printing a sparse matrix in dense format
   * \param[in] reporter is pointer to Reporter object from FaRSA
   * \param[in] name is name of Vector to print
   */
  virtual void print(const Reporter* reporter,
                     std::string name = "ans") const {};
  //@}

  /**@name get method*/
  //@{
  /**
   * get the sparsity of a matrix
   * \return is the sparsity of the given matrix
   */
  inline float sparsity() const { return sparsity_; };
  //@}
  /**@name scalar method*/
  SparseMat operator*(const double scalar);

 protected:
  /**
   * @name protected members
   */
  //@{
  unsigned int nrows_; /**number of rows*/
  unsigned int ncols_; /**number of cols*/
  std::valarray<double> val_;
  std::vector<int> rowidx_;
  std::vector<int> colidx_;
  float sparsity_;
  //@}
};  // end SparseMat

class csrMat : public SparseMat {
 public:
  csrMat(unsigned int nrows,
         unsigned int ncols,
         std::valarray<double> val,
         std::vector<int> rowidx,
         std::vector<int> colidx)
      : SparseMat(nrows, ncols, val, rowidx, colidx){};
  inline csrMat(const SparseMat& mat) : SparseMat(mat){}; /** type conversion */
  void print(const Reporter* reporter,
             std::string name = "ans") const;
  std::shared_ptr<Vector> dot(const Vector& v);
  // Vector dot(const Vector& v);
};

}  // namespace FaRSA

#endif /* __FARSASPARSEMAT_HPP__ */