/*
 * File: FaRSASparseMat.h
 * Author: Yutong Dai (rothdyt@gmail.com)
 * File Created: 2020-08-26 15:38
 * Last Modified: 2020-09-28 00:30
 * --------------------------------------------
 * Description:
 */

#ifndef __FARSASPARSEMAT_HPP__
#define __FARSASPARSEMAT_HPP__

#include <memory>
#include <string>

#include "FaRSAReporter.hpp"
#include "FaRSAVector.hpp"
namespace FaRSA {

/**
 * Forward declarations
 */
class Reporter;

class SparseMat {
 public:
  virtual void print(const Reporter* reporter,
                     std::string name) const = 0;
  virtual std::shared_ptr<Vector> dot(const Vector& v) = 0;
  virtual std::shared_ptr<SparseMat> transpose() = 0;
};

class csrMat : public SparseMat {
 public:
  csrMat(int nrows,
         int ncols,
         double* val,
         int* rowPtr,
         int* colidx)
      : nrows_(nrows),
        ncols_(ncols),
        val_(val),
        rowPtr_(rowPtr),
        colidx_(colidx){};
  csrMat(const std::string& file_path){};
  ~csrMat();
  void print(const Reporter* reporter,
             std::string name) const override;
  virtual std::shared_ptr<Vector> dot(const Vector& v) override;
  virtual std::shared_ptr<SparseMat> transpose() override;  // if I allocate new memory inside it should I mannually delete it? Where should I delete it? // introduct the notion of axes? Ax and yA might be diiferent.

 private:
  int nrows_;
  int ncols_;
  double* val_;
  int* rowPtr_;
  int* colidx_;
};

class cscMat : public SparseMat {
 public:
  cscMat(int nrows,
         int ncols,
         double* val,
         int* rowidx,
         int* colPtr)
      : nrows_(nrows),
        ncols_(ncols),
        val_(val),
        rowidx_(rowidx),
        colPtr_(colPtr){};
  cscMat(const std::string& file_path){};
  ~cscMat();
  void print(const Reporter* reporter,
             std::string name) const override;
  virtual std::shared_ptr<Vector> dot(const Vector& v) override;
  virtual std::shared_ptr<SparseMat> transpose() override;  // if I allocate new memory inside it should I mannually delete it? Where should I delete it? // introduct the notion of axes? Ax and yA might be diiferent.

 private:
  int nrows_;
  int ncols_;
  double* val_;
  int* rowidx_;
  int* colPtr_;
};

}  // namespace FaRSA

#endif /* __FARSASPARSEMAT_HPP__ */

/**
 * Abstract Sparse Matrix class
 */
// class SparseMat {
//  public:
//   /**@name Constructor*/
//   //@{
//   /**
//    * Constructor; no argument
//    */
//   SparseMat()
//       : nrows_(0),
//         ncols_(0),
//         val_({}),
//         rowidx_(0),
//         colidx_(0),
//         sparsity_(-1.0){};
//   /**
//    * Constructor with given val, rowidx, colidx
//    * \param[in] val is a pointer to an array with all nonzero elements in the SparseMat
//    * \param[in] rowidx is a pointer to an array with row indices of nonzero elements
//    * \param[in] colidx is a pointer to an array with column indices of nonzero elements
//    */
//   SparseMat(
//       unsigned int nrows,
//       unsigned int ncols,
//       std::valarray<double> val,
//       std::vector<int> rowidx,
//       std::vector<int> colidx);
//   SparseMat(const SparseMat& other)
//       : nrows_(other.nrows_),
//         ncols_(other.ncols_),
//         val_(other.val_),
//         rowidx_(other.rowidx_),
//         colidx_(other.colidx_),
//         sparsity_(other.sparsity_){};
//   //@}
//   /**@name Destructor*/
//   //@{
//   /** Destructor; val_, rowidx_, colidx_ delete*/
//   ~SparseMat();
//   //@}
//   //*@name Print  method*/
//   //@{
//   /**
//    * virtual method for printing a sparse matrix in dense format
//    * \param[in] reporter is pointer to Reporter object from FaRSA
//    * \param[in] name is name of Vector to print
//    */
//   virtual void print(const Reporter* reporter,
//                      std::string name = "ans") const {};
//   //@}

//   /**@name get method*/
//   //@{
//   /**
//    * get the sparsity of a matrix
//    * \return is the sparsity of the given matrix
//    */
//   inline float sparsity() const { return sparsity_; };
//   //@}
//   /**@name scalar method*/
//   SparseMat operator*(const double scalar);

//  protected:
//   /**
//    * @name protected members
//    */
//   //@{
//   unsigned int nrows_; /**number of rows*/
//   unsigned int ncols_; /**number of cols*/
//   std::valarray<double> val_;
//   std::vector<int> rowidx_;
//   std::vector<int> colidx_;
//   float sparsity_;
//   //@}
// };  // end SparseMat

// class csrMat : public SparseMat {
//  public:
//   csrMat(unsigned int nrows,
//          unsigned int ncols,
//          std::valarray<double> val,
//          std::vector<int> rowidx,
//          std::vector<int> colidx)
//       : SparseMat(nrows, ncols, val, rowidx, colidx){};

//   inline csrMat(const SparseMat& mat) : SparseMat(mat){}; /** type conversion for scalar multiplication(*) */

//   void print(const Reporter* reporter,
//              std::string name = "ans") const;

//   std::shared_ptr<Vector> dot(const Vector& v);
//   // Vector dot(const Vector& v);
// };