// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#include "FaRSAMatrix.hpp"

#include <stdio.h>
#include <stdlib.h>

#include "FaRSABLASLAPACK.hpp"
#include "FaRSADeclarations.hpp"

namespace FaRSA
{
// Destructor
Matrix::~Matrix()
{
    // Delete arrays
    if (column_indices_ != nullptr)
    {
        delete[] column_indices_;
        column_indices_ = nullptr;
    }  // end if
    if (row_indices_ != nullptr)
    {
        delete[] row_indices_;
        row_indices_ = nullptr;
    }  // end if
    if (values_ != nullptr)
    {
        delete[] values_;
        values_ = nullptr;
    }  // end if

}  // end destructor

// Matrix-vector product
void Matrix::matrixVectorProduct(const Vector& vector, Vector& product)
{
    // Asserts
    ASSERT_EXCEPTION(number_of_columns_ == vector.length(),
                     FARSA_MATRIX_ASSERT_EXCEPTION,
                     "Matrix assert failed.  Vector has incorrect length.");
    ASSERT_EXCEPTION(number_of_rows_ == product.length(),
                     FARSA_MATRIX_ASSERT_EXCEPTION,
                     "Matrix assert failed.  Product has incorrect length.");

    // Zero-out product
    product.scale(0.0);

    // Compute matrix-vector product, routine depending on sparse format
    if (sparse_format_ == M_COORDINATE_LIST)
    {
        for (int i = 0; i < number_of_nonzeros_; i++)
        {
            product.valuesModifiable()[row_indices_[i]] +=
                values_[i] * vector.values()[column_indices_[i]];
        }
    }
    else if (sparse_format_ == M_COMPRESSED_SPARSE_ROW)
    {
        THROW_EXCEPTION(FARSA_MATRIX_EXCEPTION,
                        "Compressed sparse row not implemented yet!");
    }  // end if
    else if (sparse_format_ == M_COMPRESSED_SPARSE_COLUMN)
    {
        THROW_EXCEPTION(FARSA_MATRIX_EXCEPTION,
                        "Compressed sparse column not implemented yet!");
    }
    else
    {
        THROW_EXCEPTION(FARSA_MATRIX_EXCEPTION, "Sparse format type error.");
    }

}  // end matrixVectorProduct

// Matrix-transpose-vector product
void Matrix::matrixTransposeVectorProduct(const Vector& vector, Vector& product)
{
    // Asserts
    ASSERT_EXCEPTION(number_of_rows_ == vector.length(),
                     FARSA_MATRIX_ASSERT_EXCEPTION,
                     "Matrix assert failed.  Vector has incorrect length.");
    ASSERT_EXCEPTION(number_of_columns_ == product.length(),
                     FARSA_MATRIX_ASSERT_EXCEPTION,
                     "Matrix assert failed.  Product has incorrect length.");

    // Zero-out product
    product.scale(0.0);

    // Compute matrix-vector product, routine depending on sparse format
    if (sparse_format_ == M_COORDINATE_LIST)
    {
        for (int i = 0; i < number_of_nonzeros_; i++)
        {
            product.valuesModifiable()[column_indices_[i]] +=
                values_[i] * vector.values()[row_indices_[i]];
        }
    }
    else if (sparse_format_ == M_COMPRESSED_SPARSE_ROW)
    {
        THROW_EXCEPTION(FARSA_MATRIX_EXCEPTION,
                        "Compressed sparse row not implemented yet!");
    }  // end if
    else if (sparse_format_ == M_COMPRESSED_SPARSE_COLUMN)
    {
        THROW_EXCEPTION(FARSA_MATRIX_EXCEPTION,
                        "Compressed sparse column not implemented yet!");
    }
    else
    {
        THROW_EXCEPTION(FARSA_MATRIX_EXCEPTION, "Sparse format type error.");
    }

}  // end matrixVectorProduct

// Set from file
void Matrix::setFromFile(char* file_name, SparseFormatType sparse_format)
{
    // Set sparse format
    sparse_format_ = sparse_format;

    // Open file
    FILE* f_in = fopen(file_name, "r");

    // Check for failed opening
    if (f_in == NULL)
    {
        THROW_EXCEPTION(FARSA_MATRIX_EXCEPTION, "Failed to open input file.");
    }

    // Read number of rows and columns (assumed first two entries in file)
    int scan_value = fscanf(f_in, "%d %d %d", &number_of_rows_,
                            &number_of_columns_, &number_of_nonzeros_);
    if (scan_value == 0 || scan_value == EOF)
    {
        THROW_EXCEPTION(FARSA_MATRIX_EXCEPTION,
                        "Number of rows and columns not read.");
    }

    // Allocate memory
    column_indices_ = new int[number_of_nonzeros_];
    row_indices_ = new int[number_of_nonzeros_];
    values_ = new double[number_of_nonzeros_];

    // Read file (assumes (row, column, value) format)
    int counter = 0;
    while (true)
    {
        int scan_value = fscanf(f_in, "%d %d %lf", &row_indices_[counter],
                                &column_indices_[counter], &values_[counter]);
        if (scan_value == 0 || scan_value == EOF ||
            counter >= number_of_nonzeros_ - 1)
        {
            break;
        }
        if (row_indices_[counter] < 0 ||
            row_indices_[counter] >= number_of_rows_)
        {
            THROW_EXCEPTION(FARSA_MATRIX_EXCEPTION, "Invalid row index read.");
        }
        if (column_indices_[counter] < 0 ||
            column_indices_[counter] >= number_of_columns_)
        {
            THROW_EXCEPTION(FARSA_MATRIX_EXCEPTION,
                            "Invalid column index read.");
        }
        counter++;
    }  // end while

    // Check if all values have been read
    if (counter < number_of_nonzeros_ - 1)
    {
        THROW_EXCEPTION(FARSA_MATRIX_EXCEPTION,
                        "Not all matrix elements have been read.");
    }

    // Close file
    fclose(f_in);

}  // end setFromFile

// Print
void Matrix::print(const Reporter* reporter, std::string name) const
{
    // Print elements of matrix
    reporter->printf(R_SOLVER, R_BASIC, "Matrix:\n");
    reporter->printf(R_SUBSOLVER, R_BASIC, "Matrix:\n");
    for (int i = 0; i < number_of_nonzeros_; i++)
    {
        reporter->printf(R_SOLVER, R_BASIC, "%s row_index(%8d)=%8d\n",
                         name.c_str(), i, row_indices_[i]);
        reporter->printf(R_SUBSOLVER, R_BASIC, "%s row_index(%8d)=%8d\n",
                         name.c_str(), i, row_indices_[i]);
    }  // end for
    for (int i = 0; i < number_of_nonzeros_; i++)
    {
        reporter->printf(R_SOLVER, R_BASIC, "%s column_index(%8d)=%8d\n",
                         name.c_str(), i, column_indices_[i]);
        reporter->printf(R_SUBSOLVER, R_BASIC, "%s column_index(%8d)=%8d\n",
                         name.c_str(), i, column_indices_[i]);
    }  // end for
    for (int i = 0; i < number_of_nonzeros_; i++)
    {
        reporter->printf(R_SOLVER, R_BASIC, "%s value(%8d)=%+23.16e\n",
                         name.c_str(), i, values_[i]);
        reporter->printf(R_SUBSOLVER, R_BASIC, "%s value(%8d)=%+23.16e\n",
                         name.c_str(), i, values_[i]);
    }  // end for

}  // end print

// col
void Matrix::col(const std::vector<int>& col_indicies, Matrix& submatrix)
{
    int counter = 0;
    if (sparse_format_ == M_COORDINATE_LIST)
    {
        for (int i = 0; i < col_indicies.size(); i++)
        {
            for (int j = 0; j < number_of_nonzeros_; j++)
            {
                // find desired colidx
                if (column_indices_[j] == col_indicies[i])
                {
                    // store the row index, column index, and values
                    submatrix.rowIndiciesModifiable()[counter] =
                        row_indices_[j];
                    submatrix.columnIndiciesModifiable()[counter] = i;
                    submatrix.valuesModifiable()[counter] = values_[j];
                    counter += 1;
                }
            }
        }  // end for
        submatrix.setSparseFormat(sparse_format_);
        submatrix.setNumberOfNonzeros(counter);
    }
    else if (sparse_format_ == M_COMPRESSED_SPARSE_ROW)
    {
        THROW_EXCEPTION(FARSA_MATRIX_EXCEPTION,
                        "Compressed sparse row not implemented yet!");
    }  // end if
    else if (sparse_format_ == M_COMPRESSED_SPARSE_COLUMN)
    {
        THROW_EXCEPTION(FARSA_MATRIX_EXCEPTION,
                        "Compressed sparse column not implemented yet!");
    }
    else
    {
        THROW_EXCEPTION(FARSA_MATRIX_EXCEPTION, "Sparse format type error.");
    }

}  // end col

}  // namespace FaRSA
