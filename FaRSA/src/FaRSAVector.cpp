// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#include <cassert>
#include <cmath>

#include "FaRSAVector.hpp"

namespace FaRSA
{

// Constructor with given length; values initialized to zero
Vector::Vector(int length)
  : length_(length),
    max_computed_(true),
    min_computed_(true),
    norm1_computed_(true),
    norm2_computed_(true),
    normInf_computed_(true),
    max_value_(0.0),
    min_value_(0.0),
    norm1_value_(0.0),
    norm2_value_(0.0),
    normInf_value_(0.0)
{

  // Allocate array
  values_ = new double[length];

  // Set values to zero
  for (int i = 0; i < length; i++) {
    values_[i] = 0.0;
  }

} // end constructor

// Constructor with given length; values initialized to given value
Vector::Vector(int length,
               double value)
  : length_(length),
    max_computed_(true),
    min_computed_(true),
    norm1_computed_(true),
    norm2_computed_(true),
    normInf_computed_(true),
    max_value_(value),
    min_value_(value)
{

  // Allocate array
  values_ = new double[length];

  // Set values to zero
  for (int i = 0; i < length; i++) {
    values_[i] = value;
  }

  // Compute norms
  norm1_value_ = (double)length * fabs(value);
  norm2_value_ = sqrt((double)length * pow(value,2.0));
  normInf_value_ = fabs(value);

} // end constructor

// Destructor; values array deleted
Vector::~Vector()
{

  // Delete array
  if (values_ != nullptr) {
    delete[] values_;
    values_ = nullptr;
  } // end if

} // end destructor

// Print array with given name
void Vector::print(const Reporter* reporter,
                   std::string name) const
{

  // Print elements
  for (int i = 0; i < length_; i++) {
    reporter->printf(R_SOLVER, R_BASIC, "%s[%6d]=%+23.16e\n", name.c_str(), i, values_[i]);
  }

} // end print

// Make new Vector as a copy
std::shared_ptr<Vector> Vector::makeNewCopy() const
{

  // Create new vector
  std::shared_ptr<Vector> vector(new Vector(length_));

  // Copy elements
  vector->copy(*this);

  // Return
  return vector;

} // end makeNewCopy

// Make new Vector by adding "scalar1" times this Vector to "scalar2" times other_vector
std::shared_ptr<Vector> Vector::makeNewLinearCombination(double scalar1,
                                                         double scalar2,
                                                         const Vector& other_vector) const
{

  // Create new vector
  std::shared_ptr<Vector> vector(new Vector(length_));

  // Copy + add elements
  vector->linearCombination(scalar1, *this, scalar2, other_vector);

  // Return
  return vector;

} // end makeNewLinearCombination

// Set length and initialize values to zero
void Vector::setLength(int length)
{

  // Store length
  length_ = length;

  // Delete previous array, if exists
  if (values_ != nullptr) {
    delete[] values_;
    values_ = nullptr;
  } // end if

  // Allocate array
  values_ = new double[length];

  // Set values to zero
  for (int i = 0; i < length_; i++) {
    values_[i] = 0.0;
  }

  // Compute scalar values
  max_computed_ = true;
  min_computed_ = true;
  norm1_computed_ = true;
  norm2_computed_ = true;
  normInf_computed_ = true;
  max_value_ = 0.0;
  min_value_ = 0.0;
  norm1_value_ = 0.0;
  norm2_value_ = 0.0;
  normInf_value_ = 0.0;

} // end setLength

// Set element with given index to given value
void Vector::set(int index,
                 double value)
{

  // Asserts
  assert(index >= 0);
  assert(index < length_);

  // Set value
  values_[index] = value;

  // Reset scalar value bools
  max_computed_ = false;
  min_computed_ = false;
  norm1_computed_ = false;
  norm2_computed_ = false;
  normInf_computed_ = false;

} // end set

// Copy elements of other_vector
void Vector::copy(const Vector& other_vector)
{

  // Assert
  assert(length_ == other_vector.length());

  // Copy elements
  for (int i = 0; i < length_; i++) {
    values_[i] = other_vector.values()[i];
  }

  // Reset scalar value bools
  max_computed_ = false;
  min_computed_ = false;
  norm1_computed_ = false;
  norm2_computed_ = false;
  normInf_computed_ = false;

} // end copy

// Copy elements of double array
void Vector::copyArray(double* array)
{

  // Copy elements
  for (int i = 0; i < length_; i++) {
    values_[i] = array[i];
  }

  // Reset scalar value bools
  max_computed_ = false;
  min_computed_ = false;
  norm1_computed_ = false;
  norm2_computed_ = false;
  normInf_computed_ = false;

} // end copyArray

// Scale elements by given scalar
void Vector::scale(double scalar)
{

  // Check for zero scalar
  if (scalar == 0.0) {

    // Scale elements
    for (int i = 0; i < length_; i++) {
      values_[i] = 0.0;
    }

  } // end if
  else if (scalar != 1.0) {

    // Scale elements
    for (int i = 0; i < length_; i++) {
      values_[i] = scalar*values_[i];
    }

  } // end else

  // Compute scalar values
  if (max_computed_) {
    if (scalar >= 0.0) {
      max_value_ = scalar*max_value_;
    }
    else {
      if (min_computed_) {
        max_value_ = scalar*min_value_;
      }
      else {
        max_computed_ = false;
      }
    }
  } // end if
  if (min_computed_) {
    if (scalar >= 0.0) {
      min_value_ = scalar*min_value_;
    }
    else {
      if (max_computed_) {
        min_value_ = scalar*max_value_;
      }
      else {
        min_computed_ = false;
      }
    }
  } // end if
  if (norm1_computed_) {
    norm1_value_ = fabs(scalar)*norm1_value_;
  }
  if (norm2_computed_) {
    norm2_value_ = fabs(scalar)*norm2_value_;
  }
  if (normInf_computed_) {
    normInf_value_ = fabs(scalar)*normInf_value_;
  }

} // end scale

// Add to this Vector "scalar" times other_vector
void Vector::addScaledVector(double scalar,
                             const Vector& other_vector)
{

  // Assert
  assert(length_ == other_vector.length());

  // Add scaled vector
  for (int i = 0; i < length_; i++) {
    values_[i] += scalar*other_vector.values()[i];
  }

  // Reset scalar value bools
  max_computed_ = false;
  min_computed_ = false;
  norm1_computed_ = false;
  norm2_computed_ = false;
  normInf_computed_ = false;

} // end addScaledVector

// Set values as linear combination (scalar1*vector1 + scalar2*vector2)
void Vector::linearCombination(double scalar1,
                               const Vector& vector1,
                               double scalar2,
                               const Vector& vector2)
{

  // Asserts
  assert(length_ == vector1.length());
  assert(length_ == vector2.length());

  // Check for nonzero scalars
  if (scalar1 != 0.0 && scalar2 != 0.0) {

    // Set elements
    for (int i = 0; i < length_; i++) {
      values_[i] = scalar1*vector1.values()[i] + scalar2*vector2.values()[i];
    }

  } // end if
  else if (scalar1 != 0.0) {

    // Set elements
    for (int i = 0; i < length_; i++) {
      values_[i] = scalar1*vector1.values()[i];
    }

  } // end else if
  else {

    // Set elements
    for (int i = 0; i < length_; i++) {
      values_[i] = scalar2*vector2.values()[i];
    }

  } // end else

  // Reset scalar value bools
  if (scalar1 == 0.0 && scalar2 == 0.0) {
    max_computed_ = true;
    min_computed_ = true;
    norm1_computed_ = true;
    norm2_computed_ = true;
    normInf_computed_ = true;
    max_value_ = 0.0;
    min_value_ = 0.0;
    norm1_value_ = 0.0;
    norm2_value_ = 0.0;
    normInf_value_ = 0.0;
  }
  else {
    max_computed_ = false;
    min_computed_ = false;
    norm1_computed_ = false;
    norm2_computed_ = false;
    normInf_computed_ = false;
  }

} // end linearCombination

// Inner product with other_vector
double Vector::innerProduct(const Vector& other_vector) const
{

  // Assert
  assert(length_ == other_vector.length());

  // Compute inner product
  double inner_product = 0.0;
  for (int i = 0; i < length_; i++) {
    inner_product += values_[i] * other_vector.values()[i];
  }

  // Return inner product
  return inner_product;

} // end innerProduct

// Maximum element
double Vector::max()
{

  // Check if computed
  if (!max_computed_) {

    // Initialize maximum
    max_value_ = values_[0];

    // Determine maximum
    for (int i = 1; i < length_; i++) {
      max_value_ = fmax(max_value_, values_[i]);
    }

    // Set to computed
    max_computed_ = true;

  } // end if

  // Return maximum
  return max_value_;

} // end max

// Minimum element
double Vector::min()
{

  // Check if computed
  if (!min_computed_) {

    // Initialize minimum
    min_value_ = values_[0];

    // Determine minimum
    for (int i = 1; i < length_; i++) {
      min_value_ = fmin(min_value_, values_[i]);
    }

    // Set to computed
    min_computed_ = true;

  } // end if

  // Return minimum
  return min_value_;

} // end max

// 1-norm
double Vector::norm1()
{

  // Check if computed
  if (!norm1_computed_) {

    // Initialize 1-norm
    norm1_value_ = fabs(values_[0]);

    // Determine 1-norm
    for (int i = 1; i < length_; i++) {
      norm1_value_ += fabs(values_[i]);
    }

    // Set to computed
    norm1_computed_ = true;

  } // end if

  // Return 1-norm
  return norm1_value_;

} // end norm1

// 2-norm
double Vector::norm2()
{

  // Check if computed
  if (!norm2_computed_) {

    // Initialize 2-norm
    norm2_value_ = pow(values_[0],2.0);

    // Determine 2-norm
    for (int i = 1; i < length_; i++) {
      norm2_value_ += pow(values_[i],2.0);
    }
    norm2_value_ = sqrt(norm2_value_);

    // Set to computed
    norm2_computed_ = true;

  } // end if

  // Return 2-norm
  return norm2_value_;

} // end norm2

// inf-norm
double Vector::normInf()
{

  // Check if computed
  if (!normInf_computed_) {

    // Initialize inf-norm
    normInf_value_ = fabs(values_[0]);

    // Determine inf-norm
    for (int i = 1; i < length_; i++) {
      normInf_value_ = fmax(normInf_value_,fabs(values_[i]));
    }

    // Set to computed
    normInf_computed_ = true;

  } // end if

  // Return inf-norm
  return normInf_value_;

} // end normInf

} // namespace NonOpt
