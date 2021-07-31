// Copyright (C) 2020 Frank E. Curtis, Daniel P. Robinson
//
// This code is published under the ??? License.
//
// Author(s) : Frank E. Curtis, Daniel P. Robinson

#include <cmath>

#include "LogisticRegression.hpp"

// Constructor
LogisticRegression::LogisticRegression(char* features_file,
                                       char* labels_file,
                                       char* groups_file,
                                       char* initial_point_file)
{

  // Read feature data
  features_.setFromFile(features_file, M_COORDINATE_LIST);

  // Read label data
  labels_.setFromFile(labels_file);

  // Read group data
  setGroupsFromFile(groups_file);

  // Read initial point
  initial_point_.setFromFile(initial_point_file);

  // Set numbers of variables and data points
  number_of_variables_ = features_.numberOfColumns();
  number_of_data_points_ = features_.numberOfRows();

} // end constructor

// Destructor
LogisticRegression::~LogisticRegression(){}

// Initial point
bool LogisticRegression::initialPoint(double* x)
{

  // Set initial point
  for (int i = 0; i < number_of_variables_; i++) {
    x[i] = initial_point_.values()[i];
  }

  // Return
  return true;

} // end initialPoint

// Objective value
bool LogisticRegression::evaluateObjective(const double* x,
                                           double& f)
{

  // Evaluate function
  f = 0.0;

  // Return
  return true;

} // end evaluateObjective

// Gradient value
bool LogisticRegression::evaluateGradient(const double* x,
                                          double* g)
{

  // Evaluate gradient
  int count = 0;
  for (int i = 0; i < groups_.size(); i++) {
    for (int j = 0; j < groups_[i].size(); j++) {
      g[count] = 0.0;
      count++;
    } // end for
  } // end for

  // Return
  return true;

} // end evaluateGradient

// Hessian-vector product
bool LogisticRegression::evaluateHessianVectorProduct(const double* x,
                                                      const std::vector<int> groups,
                                                      const double* v,
                                                      double* Hv)
{

  // Evaluate product
  int count = 0;
  for (int i = 0; i < groups.size(); i++) {
    for (int j = 0; j < groups_.at(groups.at(i)).size(); j++) {
      Hv[count] = 0.0;
      count++;
    } // end for
  } // end for

  // Return
  return true;

} // end evaluateHessianVectorProduct

// Finalize solution
bool LogisticRegression::finalizeSolution(const double* x,
                                          double f,
                                          const double* g)
{
  return true;
}

// Set groups from file
void LogisticRegression::setGroupsFromFile(char* groups_file)
{

}
