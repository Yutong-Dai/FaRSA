// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#ifndef __FARSAMATRIX_HPP__
#define __FARSAMATRIX_HPP__

#include <vector>

#include "FaRSAReporter.hpp"
#include "FaRSAVector.hpp"

namespace FaRSA
{
/**
 * Forward declarations
 */
class Reporter;
class Vector;

/**
 * Matrix class
 */
class Matrix
{
   public:
    /** @name Constructors */
    //@{
    /**
     * Constructor
     */
    Matrix()
        : number_of_columns_(0),
          number_of_nonzeros_(0),
          number_of_rows_(0),
          column_indices_(nullptr),
          row_indices_(nullptr),
          values_(nullptr){};
    /**
     * Constructor
     */
    Matrix(int number_of_rows, int number_of_columns, int number_of_nonzeros)
        : number_of_columns_(number_of_columns),
          number_of_nonzeros_(number_of_nonzeros),
          number_of_rows_(number_of_rows),
          column_indices_(new int[number_of_nonzeros]),
          row_indices_(new int[number_of_nonzeros]),
          values_(new double[number_of_nonzeros]){};
    //@}

    /** @name Destructor */
    //@{
    /**
     * Destructor
     */
    ~Matrix();
    //@}

    /** @name Print methods */
    //@{
    /**
     * Print matrix
     * \param[in] reporter is pointer to Reporter object from FaRSA
     * \param[in] name is name of Matrix to print
     */
    void print(const Reporter* reporter, std::string name) const;
    //@}

    /** @name Get methods */
    //@{
    /**
     * Get product of matrix with vector
     * \param[in] vector is reference to a Vector
     * \param[out] product is Vector to store product values
     */
    void matrixVectorProduct(const Vector& vector, Vector& product);
    /**
     * Get product of matrix transpose with vector
     * \param[in] vector is reference to a Vector
     * \param[out] product is Vector to store product values
     */
    void matrixTransposeVectorProduct(const Vector& vector, Vector& product);
    /**
     * Get number of columns
     * \return number of columns of the matrix
     */
    inline int const numberOfColumns() const { return number_of_columns_; };
    /**
     * Get number of nonzeros
     * \return number of nonzeros in the matrix
     */
    inline int const numberOfNonzeros() const { return number_of_nonzeros_; };
    /**
     * Get number of rows
     * \return number of rows of the matrix
     */
    inline int const numberOfRows() const { return number_of_rows_; };
    /**
     * @brief return values_ for modification
     *
     * \return double*
     */
    inline double* valuesModifiable() { return values_; };
    /**
     * @brief return row_indices_ for modification
     *
     * \return int*
     */
    inline int* rowIndiciesModifiable() { return row_indices_; };
    /**
     * @brief return column_indices_ for modification
     *
     * \return int*
     */
    inline int* columnIndiciesModifiable() { return column_indices_; }

    inline SparseFormatType sparseFormatType() { return sparse_format_; }
    //@}

    /** @name Set methods */
    //@{
    /**
     * Set matrix from file, compressed sparse column format
     */
    void setFromFile(char* file_name, SparseFormatType sparse_format);
    /**
     * @brief Set the sparse_format_
     *
     * \param sparse_format
     */
    inline void setSparseFormat(SparseFormatType sparse_format)
    {
        sparse_format_ = sparse_format;
    };
    /**
     * @brief Set the number_of_nonzeros_
     *
     * \param number_of_nonzeros
     */
    inline void setNumberOfNonzeros(int number_of_nonzeros)
    {
        number_of_nonzeros_ = number_of_nonzeros;
    };

    //@}

    /** @name Subsetting methods*/
    //@{

    void col(const std::vector<int>& col_indicies, Matrix& submatrix);
    //@}
   private:
    /** @name Default compiler generated methods
     * (Hidden to avoid implicit creation/calling.)
     */
    //@{
    /**
     * Copy constructor
     */
    Matrix(const Matrix&);
    /**
     * Overloaded equals operator
     */
    void operator=(const Matrix&);
    //@}

    /** @name Private members */
    //@{
    int              number_of_columns_;  /**< Number of rows of matrix */
    int              number_of_nonzeros_; /**< Number of nonzeros in matrix */
    int              number_of_rows_;     /**< Number of nonzeros in matrix */
    int*             column_indices_;     /**< Column indices */
    int*             row_indices_;        /**< Row indices */
    double*          values_;             /**< Nonzero values in matrix */
    SparseFormatType sparse_format_;      /**< Sparse format type */
                                          //@}

};  // end Matrix

}  // namespace FaRSA

#endif /* __FARSAMATRIX_HPP__ */
